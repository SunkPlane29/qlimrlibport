## About This Project

This project is a refactored and adapted version of the original C++ code MUSES/QLIMR, made by Carlos Conde. The original source can be found here: https://gitlab.com/nsf-muses/ns-qlimr/qlimr.

The primary purpose of this adaptation was to convert the original cluster/container-based application into a portable shared library (`.so`, `.dll`) with a C-style API. This makes the core logic accessible to other programming languages (like Julia, Python, etc.) and simplifies its use by replacing file-based inputs with direct function calls.

**Key Changes:**
* Refactored the core logic into a shared library.
* Created a C-style `extern "C"` API for cross-language compatibility.
* Removed dependencies on specific file structures for input/output.
* Simplified the interface to a set of direct function calls.

## License

This project is licensed under the **GNU General Public License v3.0** (or the same version as the original).


This library is an adaptation of the original work MUSES/QLIMR by Carlos Conde, which is also licensed under the GPL. As a derivative work, this project is bound by the terms of the original license.
