#pragma once
#define LOGGING_LEVEL_RELEASE 0
#define LOGGING_LEVEL_DEBUG 1
#define LOGGING_LEVEL_TRACE 2

#ifndef LOGGING_LEVEL
#define LOGGING_LEVEL LOGGING_LEVEL_RELEASE
#endif 

#if LOGGING_LEVEL > LOGGING_LEVEL_RELEASE
#include <string>

class PrintLog {
public:
	PrintLog& operator()();
	PrintLog& operator()(const char*);
	PrintLog& operator()(const std::string&);
	PrintLog& operator()(float);
	PrintLog& operator()(long long);
	PrintLog& operator()(int);
	PrintLog& operator()(unsigned);
	PrintLog& operator()(size_t);

	static PrintLog& instance() {
		static PrintLog theInstance;
		return theInstance;
	}
	void init(const std::string& experimentName);
	void init_stdout();
private:
	PrintLog() = default;

	void printlog_part(const std::string& line_part);
	void printlog_end_of_line();
};

inline PrintLog& printlog() {
	return PrintLog::instance()();
}
template<typename T>
PrintLog& printlog(T&& value) {
	return PrintLog::instance()(value);
}

#if FULL_STAT_LOGGING >= LOGGING_LEVEL_DEBUG
inline PrintLog& printlog_debug() { 
	return PrintLog::instance()();
}
template<typename T>
PrintLog& printlog_debug(T&& value) {
	return PrintLog::instance()(value);
}
#else
#define printlog_debug(...) DontPrintLog()()
#endif

#if FULL_STAT_LOGGING >= LOGGING_LEVEL_TRACE
inline PrintLog& printlog_trace() {
	return PrintLog::instance()();
}
template<typename T>
PrintLog& printlog_trace(T&& value) {
	return PrintLog::instance()(value);
}
#else
#define printlog_trace(...) DontPrintLog()()
#endif


inline void init_logger(const std::string& experimentName) {
	PrintLog::instance().init(experimentName);
}
inline void init_logger() {
	PrintLog::instance().init_stdout();
}
#else
struct DontPrintLog {
	DontPrintLog& operator()(...) { return *this; }
};

#define init_logger(...)
#define printlog(...) DontPrintLog()()
#define printlog_debug printlog
#define printlog_trace printlog
#endif // FULL_STAT_LOGGING