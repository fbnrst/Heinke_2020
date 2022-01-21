#include <stdio.h>
#include <list>
#include <gsl/gsl_interp.h>
#include <limits>
#include <random>

enum cell_type { t_2n, t_2x2n, t_4n };
const unsigned int NUMCELLS =3;
enum event { die, live };
class parameter_T{
   public:
   parameter_T(){}
  ~parameter_T(){
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
class parameter{
  public:
  parameter(void):
    bernoulli_d(NULL),
    random_generator(NULL),
    value(0){}
  parameter(double val):
    bernoulli_d(NULL),
    random_generator(NULL),
    value(val){}
  ~parameter(){
    delete bernoulli_d;
  }
  double operator() (void){
    return (*bernoulli_d)(*random_generator)*value;
  }
  
  void set_rand(std::mt19937 *r_gen,double dt){
    this->dt = dt;
    random_generator = r_gen;
    delete bernoulli_d;
    bernoulli_d = new std::bernoulli_distribution(value*dt);
  }
  void update_value(double val){
    value = val;
    delete bernoulli_d;
    bernoulli_d = new std::bernoulli_distribution(value*dt);
  }
  std::bernoulli_distribution *bernoulli_d;
  std::mt19937 *random_generator;
  double value;
  double dt;
};

class parameters {
 public:
  parameters(){}
  ~parameters(){
  }
  void set_random(std::mt19937 *r_gen,double dt){
    gen = r_gen;
    this->dt = dt;
    d2n.set_rand(r_gen,dt);
    d4n.set_rand(r_gen,dt);
    d2x2n.set_rand(r_gen,dt);
    g2x2nb2n.set_rand(r_gen,dt);
    k2nb2x2n.set_rand(r_gen,dt);
    k2nb4n.set_rand(r_gen,dt);

    g2nb2n.set_rand(r_gen,dt);
    k2x2nb4n.set_rand(r_gen,dt);
    k4nb4n.set_rand(r_gen,dt);
    k4nb2x2n.set_rand(r_gen,dt);
    k2x2nb2x2n.set_rand(r_gen,dt);
  }
  void get(double now){
    c_time = now;
    g2nb2n.update_value(g2nb2n_t.get(now));
    k2x2nb4n.update_value(k2x2nb4n_t.get(now));
    k4nb4n.update_value(k4nb4n_t.get(now));
    k4nb2x2n.update_value(k4nb2x2n_t.get(now));
    k2x2nb2x2n.update_value(k2x2nb2x2n_t.get(now));
    }
  std::mt19937 *gen;
  double dt;

  parameter_T g2nb2n_t;
  parameter_T k2x2nb4n_t;
  parameter_T k4nb4n_t;
  parameter_T k4nb2x2n_t;
  parameter_T k2x2nb2x2n_t;
  double max_age;
  double c_time;
  parameter d2n;
  parameter d4n;
  parameter d2x2n;
  parameter g2x2nb2n;
  parameter k2nb2x2n;
  parameter k2nb4n;

  parameter g2nb2n;
  parameter k2x2nb4n;
  parameter k4nb4n;
  parameter k4nb2x2n;
  parameter k2x2nb2x2n;
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

class C2x2N: public cell_base {
  public:
    C2x2N(parameters &A, double btime) : cell_base(A,btime,t_2x2n) {}
    ~C2x2N() {}
    void events(std::list<cell_base*>::iterator &c_it,std::list<cell_base*> &current_cells, std::list<cell_base*> &old_cells);
};

class C4N: public cell_base {
  public:
    C4N(parameters &A, double btime) : cell_base(A,btime,t_4n) {}
    ~C4N() {}
    void events(std::list<cell_base*>::iterator &c_it,std::list<cell_base*> &current_cells, std::list<cell_base*> &old_cells);
};


class sim{
  public:
  std::vector<std::vector<std::vector<int>>> time_age;
  void run(parameters &p, double delta_t, double delta_t_out, double maxT, unsigned int seed, unsigned int num_2n, unsigned int num_2x2n, unsigned int num_4n);
  std::vector<std::vector<std::vector<int>>> get_time_age(){return time_age;};
};   

void check_cell(double OODT, cell_base *cell, std::vector<std::vector<std::vector<int>>>& time_age);
