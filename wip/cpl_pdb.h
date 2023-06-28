/**
 * @file cpl_pdb.h
 *
 * @brief This file contains the API for reading protein data bank files.
 * */

#pragma once

#ifndef CPL_PDB_H
#define CPL_PDB_H

#if !defined(CPL_PDB_FUNC) && defined(CPL_PDB_STATIC)
#define CPL_PDB_FUNC static
#endif

#if !defined(CPL_PDB_FUNC) && defined(__cplusplus)
#define CPL_PDB_FUNC extern "C"
#endif

#ifndef CPL_PDB_FUNC
#define CPL_PDB_FUNC
#endif

/**
 * @defgroup cpl_pdb_api PDB API
 *
 * @brief An API for reading protein data bank files.
 * */

/**
 * @brief This structure contains pointers to the callback functions that will forward the data to the client code.
 *
 * @note The functions in this structure may be called more than once.
 *
 * @ingroup cpl_pdb_api
 * */
struct cpl_pdb_reader
{
  /** @brief Called when a portion of the header is read. */
  void (*on_header)(void* user_data, const char* header, int header_size);

  /** @brief Called when an error occurs while reading the file. */
  void (*on_error)(void* user_data, const char* msg);
};

#ifndef CPL_PDB_NO_STDIO

/**
 * @brief Reads a PDB file from the file system.
 *
 * @param path The path to the PDB file to read.
 *
 * @param user_data An optional pointer to pass to the callback functions.
 *
 * @param reader A structure containing pointers to the various callbacks to pass the data to.
 *
 * @return Zero on success, negative one on failure.
 *
 * @ingroup cpl_pdb_api
 * */
CPL_PDB_FUNC int
cpl_pdb_read(const char* path, void* user_data, const struct cpl_pdb_reader* reader);

#endif /* CPL_PDB_NO_STDIO */

CPL_PDB_FUNC int
cpl_pdb_read_from_stream(void* stream_data,
                         int (*stream_read_func)(void* stream_data_ptr, char* buf, int size),
                         void* user_data,
                         const struct cpl_pdb_reader* reader);

CPL_PDB_FUNC int
cpl_pdb_read_from_line(const char* line, int line_size, void* user_data, const struct cpl_pdb_reader* reader);

#ifdef CPL_PDB_IMPLEMENTATION

#define CPL_PDB_ERROR(msg)                                                                                             \
  do {                                                                                                                 \
    if (reader->on_error)                                                                                              \
      reader->on_error(user_data, (msg));                                                                              \
  } while (0)

#ifndef CPL_PDB_NO_STDIO
#include <stdio.h>
#endif

#ifndef CPL_PDB_NO_STDIO

static int
cpl_pdb_stdio_read(void* file_ptr, char* buf, int size)
{
  FILE* file;

  size_t read_size;

  file = file_ptr;

  if (feof(file))
    return 0;

  /* for testing purposes, max out read size to a value of 3
   * remove this after testing
   * */
  size = (size > 3) ? 3 : size;

  read_size = fread(buf, 1, (size_t)size, file);

  if (read_size != ((size_t)size)) {
    if (!feof(file))
      return -1;
  }

  return (int)read_size;
}

CPL_PDB_FUNC int
cpl_pdb_read(const char* path, void* user_data, const struct cpl_pdb_reader* reader)
{
  FILE* file;

  int result;

  file = fopen(path, "rb");
  if (!file) {
    CPL_PDB_ERROR("failed to open file");
    return -1;
  }

  result = cpl_pdb_read_from_stream(file, cpl_pdb_stdio_read, user_data, reader);

  (void)fclose(file);

  return result;
}

#endif /* CPL_PDB_NO_STDIO */

/**
 * @defgroup cpl_pdb_internals PDB Internals
 *
 * @brief The internal documentation of the PDB reader API.
 * */

#ifndef CPL_PDB_COLUMN_LIMIT
#define CPL_PDB_COLUMN_LIMIT 80
#endif

#ifndef CPL_PDB_READ_MAX
#define CPL_PDB_READ_MAX (CPL_PDB_COLUMN_LIMIT * 100)
#endif

/**
 * @brief Locates the end of the line in a buffer, if it exists.
 *
 * @return The index of the character immediately *after* the newline sequence.
 *
 * @ingroup cpl_pdb_internals
 * */
CPL_PDB_FUNC int
cpl_pdb__find_eol(const char* buffer, int buffer_size)
{
  int i;

  for (i = 0; i < buffer_size; i++) {
    if (buffer[i] == '\n')
      return i + 1;

    if (((i + 1) < buffer_size) && (buffer[i] == '\r') && (buffer[i + 1] == '\n'))
      return i + 2;
  }

  return -1;
}

CPL_PDB_FUNC int
cpl_pdb_read_from_stream(void* stream_data,
                         int (*stream_read_func)(void* stream_data_ptr, char* buf, int size),
                         void* user_data,
                         const struct cpl_pdb_reader* reader)
{
  char read_buffer[CPL_PDB_READ_MAX];
  int total_read_size;
  int read_size;
  int line_length;
  int offset;
  int offset2;
  int remaining;
  int eof;

  total_read_size = 0;

  eof = 0;

  while (!eof) {

    /* read until buffer is full or the end of file is reached. */

    while (total_read_size < CPL_PDB_READ_MAX) {

      read_size = stream_read_func(stream_data, &read_buffer[total_read_size], CPL_PDB_READ_MAX - total_read_size);

      if (read_size == 0) {
        eof = 1;
        break;
      }

      if (read_size == -1) {
        CPL_PDB_ERROR("failed to read from stream");
        return -1;
      }

      total_read_size += read_size;
    }

    /* parse all available lines in the buffer */

    offset = 0;

    while (offset < total_read_size) {

      remaining = total_read_size - offset;

      line_length = cpl_pdb__find_eol(read_buffer + offset, remaining);

      if (line_length == -1) {
        if (remaining >= CPL_PDB_COLUMN_LIMIT) {
          CPL_PDB_ERROR("column limit exceeded");
          return -1;
        }
        break;
      }

      if (cpl_pdb_read_from_line(read_buffer + offset, line_length, user_data, reader) == -1)
        return -1;

      offset += line_length;
    }

    /* shift whatever was not able to be parsed to the beginning of the read buffer */

    remaining = total_read_size - offset;

    for (offset2 = 0; offset2 < remaining; offset2++)
      read_buffer[offset2] = read_buffer[offset + offset2];

    total_read_size = remaining;

    if (eof && (total_read_size > 0)) {
      CPL_PDB_ERROR("last line is missing newline sequence");
      return -1;
    }
  }

  return 0;
}

CPL_PDB_FUNC int
cpl_pdb_read_from_line(const char* line, int line_size, void* user_data, const struct cpl_pdb_reader* reader)
{
  (void)line;
  (void)line_size;
  (void)user_data;
  (void)reader;
  return 0;
}

#endif /* CPL_PDB_IMPLEMENTATION */

#endif /* CPL_PDB_H */
