#pragma once

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

class LogComparer
{
public:
  // Constructor
  LogComparer(std::string ref_log, std::string new_log)
      : ref_log_file(std::move(ref_log))
      , new_log_file(std::move(new_log))
  {
  }

  // Supply any *literal string* patterns to ignore lines containing those
  // strings
  void setIgnorePatterns(const std::vector<std::string>& patterns)
  {
    ignore_patterns = patterns;
  }

  // Run the diff with filtering
  auto compareLogs() const -> bool
  {
    const std::string diff_cmd = buildDiffCommand();
    int diff_status = system(diff_cmd.c_str());
    return diff_status == 0;  // Return true if logs match, false otherwise
  }

  // Get the generated diff command (for debugging/logging)
  auto getDiffCommand() const -> std::string { return buildDiffCommand(); }

private:
  const std::string ref_log_file;
  const std::string new_log_file;
  std::vector<std::string> ignore_patterns;

  // Build the grep -F command that ignores lines containing any pattern
  auto buildGrepFilterCommand() const -> std::string
  {
    // No patterns to ignore
    if (ignore_patterns.empty()) {
      return "";
    }

    // Example:
    //   grep -Fv -e "PAT1" -e "PAT2" -e "PAT3" ...
    std::ostringstream cmd;
    cmd << "| grep -Fv";
    for (const auto& pat : ignore_patterns) {
      cmd << " -e \"" << pat << "\"";
    }
    return cmd.str();
  }

  // Build the full diff command with filtering
  auto buildDiffCommand() const -> std::string
  {
    std::string grep_filter = buildGrepFilterCommand();

    std::cout << "Comparing " << ref_log_file << " and " << new_log_file
              << std::endl;
    std::cout << "Grep filter: " << grep_filter << std::endl;

    // We'll replace memory addresses and timing info, then pipe both logs
    // through the same sed+grep filter, then compare with diff -u.
    std::ostringstream cmd;
    cmd << "bash -c 'diff -u "
        << R"(<(sed -E "s/0x[0-9a-f]+/0xADDR/g; s/[0-9]+\.[0-9]+ s/XX.XX s/g" )"
        << ref_log_file << " " << grep_filter << ") "
        << R"(<(sed -E "s/0x[0-9a-f]+/0xADDR/g; s/[0-9]+\.[0-9]+ s/XX.XX s/g" )"
        << new_log_file << " " << grep_filter << ")'";

    return cmd.str();
  }
};
