#include <gtest/gtest.h>

#include "../cpl_pdb.h"

namespace {

struct reader_impl final
{
  std::string errors;

  static void on_error(void* self_ptr, const char* msg, int)
  {
    auto* self = static_cast<reader_impl*>(self_ptr);
    self->errors += msg;
    self->errors += "\n";
  }
};

} // namespace

TEST(pdb, read)
{
  reader_impl reader_data;

  const struct cpl_pdb_reader reader
  {
    &reader_impl::on_error, nullptr,
  };

  cpl_pdb_read(PROJECT_SOURCE_DIR "/pdb/1AKI.pdb", &reader_data, &reader);

  EXPECT_EQ(reader_data.errors, "");
}
