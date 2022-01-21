#include <stdio.h>
#include <list>
#include <gsl/gsl_interp.h>
#include <limits>
#include <random>

enum cell_type { t_2n, t_4n };
enum event { die, live };
class parameter{
   public:
  ~parameter(){
      gsl_interp_free(interp);
      delete[] tp;
      delete[] r;
  }
  void set(std::list<double> timep, std::list<double> rate);
  double get(double now);
  private:
    gsl_interp *interp;
    double *tp,*r;
};


class parameters {
 public:
  parameters():
    b_lambda2(NULL),
    b_lambda4(NULL),
    b_kappa24(NULL),
    b_kappa42(NULL),
    b_delta2(NULL),
    b_delta4(NULL) {}
  ~parameters(){
    delete b_lambda2;
    delete b_lambda4;
    delete b_kappa24;
    delete b_kappa42;
    delete b_delta2;
    delete b_delta4;
  }
  void set_random(std::mt19937 *r_gen,double dt){
    delta_t = dt;
    gen = r_gen;
    delete b_kappa24;
    b_kappa24 = new std::bernoulli_distribution(kappa24*dt);
    delete b_kappa42;
    b_kappa42 = new std::bernoulli_distribution(kappa42*dt); 
    delete b_delta2;
    b_delta2 = new std::bernoulli_distribution(delta2*dt);
    delete b_delta4;
    b_delta4 = new std::bernoulli_distribution(delta4*dt);
  }
  void get(double now){
    c_time = now;
    lambda2 = lambda2_t.get(now);
    lambda4 = lambda4_t.get(now);
    delete b_lambda2;
    b_lambda2 = new std::bernoulli_distribution(lambda2*delta_t);
    delete b_lambda4;
    b_lambda4 = new std::bernoulli_distribution(lambda4*delta_t);
    }
  bool r_lambda2(){return (*b_lambda2)(*gen);}
  bool r_lambda4(){return (*b_lambda4)(*gen);}
  bool r_kappa24(){return (*b_kappa24)(*gen);}
  bool r_kappa42(){return (*b_kappa42)(*gen);}
  bool r_delta2() {return (*b_delta2) (*gen);}
  bool r_delta4() {return (*b_delta4) (*gen);}
  parameter lambda2_t;
  parameter lambda4_t;
  double max_age;
  double delta_t;
  double c_time;
  double lambda2;
  double kappa24;
  double lambda4;
  double kappa42;
  double delta2;
  double delta4;
  std::mt19937 *gen;
  std::bernoulli_distribution *b_lambda2;
  std::bernoulli_distribution *b_lambda4;
  std::bernoulli_distribution *b_kappa24;
  std::bernoulli_distribution *b_kappa42;
  std::bernoulli_distribution *b_delta2;
  std::bernoulli_distribution *b_delta4;
  private:   
};

class cell_base {
  public:
    cell_base(parameters &A, double btime, cell_type c_t) : 
        p(A),
        type(c_t),
        t_birth(btime),
        t_death(A.max_age + 100) {} 
    cell_base(parameters &&A, double btime, cell_type c_t) = delete; 
    virtual ~cell_base(){}

    virtual void events(std::list<cell_base*>::iterator &c_it,std::list<cell_base*> &current_cells, std::list<cell_base*> &old_cells){std::cout<<"virt;";return;}
    parameters &p;
    cell_type type;
    double t_birth;
    double t_death;
  private:
};


class C2N: public cell_base {
  public:
    C2N(parameters &A, double btime) : cell_base(A,btime,t_2n) {}
    ~C2N() {}
    void events(std::list<cell_base*>::iterator &c_it,std::list<cell_base*> &current_cells, std::list<cell_base*> &old_cells);
};

class C4N: public cell_base {
  public:
    C4N(parameters &A, double btime) : cell_base(A,btime,t_4n) {}
    ~C4N() {}
    void events(std::list<cell_base*>::iterator &c_it,std::list<cell_base*> &current_cells, std::list<cell_base*> &old_cells);
};

std::pair<std::vector<std::vector<int>>,std::vector<std::vector<int>>> test(parameters &p, double delta_t, double delta_t_out, double maxT, unsigned int seed, unsigned int num_2n, unsigned int num_4n);
void check_cell(double DT, double OODT, cell_base *cell, std::vector<std::vector<int>>& n2_time_age,std::vector<std::vector<int>>& n4_time_age);
