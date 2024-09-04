#pragma once

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <nlohmann/json.hpp>

// for convenience
using json = nlohmann::json;

/**
 * @class JsonComparer
 * @brief A class for comparing JSON files.
 *
 * The JsonComparer class provides functionality to compare two JSON files and
 * check for differences between them. It supports comparing objects, arrays,
 * numbers, and other value types. The comparison can be customized by setting
 * a tolerance value (percent) for numeric comparisons.
 */
class JsonComparer
{
public:
  JsonComparer(double tolerance = 0.0, int print_limit = 50)
      : m_tolerance(tolerance)
      , m_print_limit(print_limit)
  {
  }

  auto compare(const std::string& file1, const std::string& file2) -> bool
  {
    m_print_count = 0, m_output_suppressed = false;

    json json1, json2;

    if (!read_json(file1, json1) || !read_json(file2, json2)) {
      return false;
    }

    bool are_equal = compare_json(json1, json2, "");
    if (!are_equal) {
      std::cout << "Differences found:\n";
    }

    return are_equal;
  }

private:
  int m_print_limit, m_print_count;
  bool m_output_suppressed;
  double m_tolerance;

  static auto read_json(const std::string& file_path, json& json_obj) -> bool
  {
    if (!std::filesystem::exists(file_path)) {
      throw std::runtime_error("File does not exist: " + file_path);
    }

    std::ifstream ifs(file_path);
    if (!ifs.is_open()) {
      throw std::runtime_error("Error opening file: " + file_path);
    }

    try {
      ifs >> json_obj;
    } catch (const json::parse_error& e) {
      throw std::runtime_error("Parse error in file " + file_path + ": "
                               + e.what());
    }

    return true;
  }

  auto compare_json(const json& json1,
                    const json& json2,
                    const std::string& path) -> bool
  {
    if (json1.type() != json2.type()) {
      print_difference(
          "Type mismatch", json1.type_name(), json2.type_name(), path);
      return false;
    }

    switch (json1.type()) {
      case json::value_t::object:
        return compare_objects(json1, json2, path);
      case json::value_t::array:
        return compare_arrays(json1, json2, path);
      case json::value_t::number_float:
      case json::value_t::number_integer:
      case json::value_t::number_unsigned:
        return compare_numbers(json1, json2, path);
      default:
        if (json1 != json2) {
          print_difference("Value mismatch", json1, json2, path);
          return false;
        }
        break;
    }

    return true;
  }

  auto compare_objects(const json& json1,
                       const json& json2,
                       const std::string& path) -> bool
  {
    bool are_equal = true;
    for (auto it = json1.begin(); it != json1.end(); ++it) {
      std::string new_path = path.empty() ? it.key() : path + "." + it.key();
      if (json2.contains(it.key())) {
        if (!compare_json(it.value(), json2.at(it.key()), new_path)) {
          are_equal = false;
        }
      } else {
        print_difference(
            "Missing key in second JSON", it.value(), "", new_path);
        are_equal = false;
      }
    }
    for (auto it = json2.begin(); it != json2.end(); ++it) {
      std::string new_path = path.empty() ? it.key() : path + "." + it.key();
      if (!json1.contains(it.key())) {
        print_difference("Extra key in second JSON", "", it.value(), new_path);
        are_equal = false;
      }
    }
    return are_equal;
  }

  auto compare_arrays(const json& json1,
                      const json& json2,
                      const std::string& path) -> bool
  {
    bool are_equal = true;
    size_t min_size = std::min(json1.size(), json2.size());
    for (size_t i = 0; i < min_size; ++i) {
      if (!compare_json(
              json1[i], json2[i], path + "[" + std::to_string(i) + "]"))
      {
        are_equal = false;
      }
    }
    if (json1.size() > json2.size()) {
      for (size_t i = json2.size(); i < json1.size(); ++i) {
        print_difference("Extra element in first JSON",
                         json1[i],
                         "",
                         path + "[" + std::to_string(i) + "]");
        are_equal = false;
      }
    } else if (json2.size() > json1.size()) {
      for (size_t i = json1.size(); i < json2.size(); ++i) {
        print_difference("Extra element in second JSON",
                         "",
                         json2[i],
                         path + "[" + std::to_string(i) + "]");
        are_equal = false;
      }
    }
    return are_equal;
  }

  auto compare_numbers(const json& json1,
                       const json& json2,
                       const std::string& path) -> bool
  {
    double val1 = json1.get<double>();
    double val2 = json2.get<double>();
    double percent_diff =
        std::fabs(val1 - val2) / std::max(std::fabs(val1), std::fabs(val2));
    if (percent_diff > m_tolerance) {
      print_difference("Value mismatch", val1, val2, path, percent_diff * 100);
      return false;
    }
    return true;
  }

  void print_difference(const std::string& message,
                        const json& json1,
                        const json& json2,
                        const std::string& path,
                        double percent_diff = 0.0)
  {
    if (m_print_count >= m_print_limit) {
      if (!m_output_suppressed) {
        std::cout << "Suppressing further output after reaching the limit of "
                  << m_print_limit << " differences\n";
      }
      m_output_suppressed = true;
      return;
    }

    if (json1.is_null() && json2.is_null()) {
      std::cout << message << " at " << path << ": both values are null\n";
    } else if (json1.is_null() || json2.is_null()) {
      std::cout << message << " at " << path << ": " << json1 << " vs " << json2
                << "\n";
    } else {
      std::cout << message << " at " << path << ": " << json1 << " vs "
                << json2;
      if (percent_diff != 0.0) {
        std::cout << " (difference: " << percent_diff << "%)";
      }
      std::cout << "\n";
    }
    m_print_count++;
  }
};
