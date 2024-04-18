#include <gtest/gtest.h>

#include "TFile.h"

TEST(ParamFileTests, Read)
{
  const std::string input_file_path =
      std::string(TEST_DATA_DIR) + "SelEkin_id22_Mom65536_65536_eta_20_25.root";

  std::unique_ptr<TFile> param_file(
      TFile::Open(input_file_path.c_str(), "READ"));

  Long64_t bytes_read = param_file->GetBytesRead();

  EXPECT_EQ(bytes_read, 8312);
}
