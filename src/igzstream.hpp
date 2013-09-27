#ifndef __IGZSTREAM_HPP
#define __IGZSTREAM_HPP


#include <cstdlib>
#include <iostream>
#include <fstream>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>


class igzstream : public boost::iostreams::filtering_istream {
private:
  // if this is not cin, store ifstream object, and delete it when done
  std::unique_ptr<std::ifstream> p_file;

  // no copy constructor
  igzstream(const igzstream &) : std::basic_ios<char>(), boost::iostreams::filtering_istream() {}

  // no copy-assignment operator
  igzstream & operator = (const igzstream &) { return *this; }

public:
  igzstream() {}

  igzstream(const char * name) { open(name); }
  igzstream(const std::string& name) { open(name.c_str()); }

  igzstream(std::istream & is) { attach(is); }

  void open(const char * name) {
    if (strncmp(name, "-", 2) == 0) {
      attach(std::cin);
    } else {
      p_file = std::unique_ptr<std::ifstream>(new std::ifstream(name));
      if (!*p_file) {
	std::cerr << "error opening file [" << name << "]\n";
	exit(EXIT_FAILURE);
      }
      attach(*p_file);
    }
  }

  void attach(std::istream& is) {
    int c = is.peek();
    is.clear();
    if (c == 31) {
      push(boost::iostreams::gzip_decompressor());
    }
    push(is);
  }
};


#endif
