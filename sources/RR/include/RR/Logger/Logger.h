#pragma once
#include "LoggingLevel.h"

#ifndef LOGGING_LEVEL
#define LOGGING_LEVEL LOGGING_LEVEL_DEBUG
#endif // !LOGGING_LEVEL

#include <string>

namespace RR::Logger {
	struct DontPrintLog {
		template<typename ...Args>
		DontPrintLog& operator()(Args&&...) { return *this; }
	};

	class PrintLog {
	public:
		PrintLog& operator()() {
			printlog_end_of_line();
			return *this;
		}
		PrintLog& operator()(const char* line) {
			printlog_part(line);
			return *this;
		}
		template<typename T>
		PrintLog& operator()(const T& val) {
			printlog_part(std::to_string(val));
			return *this;
		}
		PrintLog& operator()(const std::string& line) {
			printlog_part(line);
			return *this;
		}

		static PrintLog& instance() {
			static PrintLog theInstance;
			return theInstance;
		}
		void init(const std::string& experiment_name);
		void init_stdout();
		bool initialized() const;
	private:
		PrintLog() = default;

		void printlog_part(const std::string& line_part);
		void printlog_end_of_line();
	};

	// functions interface:

#if LOGGING_LEVEL >= LOGGING_LEVEL_RELEASE
	inline void init_logger(const std::string& log_path) {
		PrintLog::instance().init(log_path);
	}
	inline void init_logger() {
		PrintLog::instance().init_stdout();
	}

	inline PrintLog& printlog() {
		return PrintLog::instance()();
	}
	template<typename T>
	PrintLog& printlog(T&& value) {
		return PrintLog::instance()(value);
	}
#else
	inline void init_logger(const std::string& log_path) { }
	inline void init_logger() { }

	inline DontPrintLog printlog() {
		return DontPrintLog();
	}
	template<typename T>
	DontPrintLog printlog(T&& value) {
		return DontPrintLog();
	}
#endif

#if LOGGING_LEVEL >= LOGGING_LEVEL_DEBUG
	inline PrintLog& printlog_debug() {
		return PrintLog::instance()();
	}
	template<typename T>
	PrintLog& printlog_debug(T&& value) {
		return PrintLog::instance()(value);
	}
#else
	inline DontPrintLog printlog_debug() {
		return DontPrintLog();
	}
	template<typename T>
	DontPrintLog printlog_debug(T&& value) {
		return DontPrintLog();
	}
#endif

#if LOGGING_LEVEL >= LOGGING_LEVEL_TRACE
	inline PrintLog& printlog_trace() {
		return PrintLog::instance()();
	}
	template<typename T>
	PrintLog& printlog_trace(T&& value) {
		return PrintLog::instance()(value);
	}
#else 
	inline DontPrintLog printlog_trace() {
		return DontPrintLog();
	}
	template<typename T>
	DontPrintLog printlog_trace(T&& value) {
		return DontPrintLog();
	}
#endif
}