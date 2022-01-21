%module POP3p2x2
%{
#   define SWIG_PYTHON_EXTRA_NATIVE_CONTAINERS
%}
%include <std_list.i>
%include <std_pair.i>
%include <std_vector.i>
%exception { 
    try {
        $action
    } catch (const std::exception& e) {
        SWIG_Error(SWIG_RuntimeError, e.what());
    } catch (...) {
        SWIG_exception(SWIG_RuntimeError, "unknown exception");
    }
}

%template(Listd) std::list<double> ;

%template(vecint) std::vector<int>;
%template(vecvecint) std::vector<std::vector<int> >;
%template(vecvecvecint) std::vector<std::vector<std::vector<int>  >>;

%{
  #include "POP3p2x2.h"
%}
%include "POP3p2x2.h"

