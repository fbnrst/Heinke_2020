%module Rm
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
%template(paivecvecint) std::pair<std::vector<std::vector<int> >, std::vector<std::vector<int> >>;

%{
  #include "Rm.h"
%}
%include "Rm.h"

