// Copyright (c) 2026 CERN for the benefit of the FastCaloSim project

#pragma once

#include <cstdint>

/// Result of a simulate call: success, a fatal error, or a request to simulate
/// again. A retry request carries the number of retries the caller may make,
/// e.g. FCSReturnCode::retry(5).
class FCSReturnCode
{
public:
  enum class Status : std::uint8_t
  {
    Fatal,
    Success,
    Retry
  };

  explicit constexpr FCSReturnCode(Status status) noexcept
      : m_status(status)
  {
  }

  /// A request to simulate again, at most max_retries times
  static constexpr auto retry(int max_retries) noexcept -> FCSReturnCode
  {
    FCSReturnCode code(Status::Retry);
    code.m_max_retries = max_retries;
    return code;
  }

  constexpr auto status() const noexcept -> Status { return m_status; }
  constexpr auto is_retry() const noexcept -> bool
  {
    return m_status == Status::Retry;
  }
  /// The number of retries a retry request allows, 0 for the other results
  constexpr auto max_retries() const noexcept -> int { return m_max_retries; }

  friend constexpr auto operator==(FCSReturnCode, FCSReturnCode) noexcept
      -> bool = default;

private:
  Status m_status;
  int m_max_retries = 0;
};

inline constexpr FCSReturnCode FCSSuccess {FCSReturnCode::Status::Success};
inline constexpr FCSReturnCode FCSFatal {FCSReturnCode::Status::Fatal};
