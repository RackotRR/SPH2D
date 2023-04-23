#include <RR/Logger/Logger.h>

#include <fstream>
#include <iostream>
#include <memory>

using namespace RR::Logger;

static std::unique_ptr<std::ofstream> log_file;
static std::ostream* log_stream = nullptr;

void PrintLog::init(const std::string& log_path) {
	std::string filename = log_path;
	log_file = std::make_unique<std::ofstream>(filename);
	log_stream = log_file.get();
}
void PrintLog::init_stdout() {
	log_stream = &std::cout;
}
void PrintLog::printlog_part(const std::string& line_part) {
	if (initialized()) *log_stream << line_part;
}
void PrintLog::printlog_end_of_line() {
	if (initialized()) *log_stream << std::endl;
}

bool PrintLog::initialized() const {
	return log_file || log_stream;
}