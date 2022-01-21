#include <stdio.h>
#include <iostream>
#include <list>
#include "POP3p2x2.h"
#include <tuple>
#include <stdexcept>
#include <gsl/gsl_interp.h>
#define SLOW
void parameter_T::set(std::list<double> timep,std::list<double> rate){
    if (timep.size() != rate.size()){
        throw std::invalid_argument(  "time and rate not same length" );   
        return;
    }
    unsigned int size = rate.size();
    tp = new double[size];
    r = new double[size];
    unsigned int pos=0;
    for (auto elm : timep){
        tp[pos] = elm;
        pos++;
    }
    pos=0;
    for (auto elm : rate){
        r[pos] = elm;
        pos++;
    }
    interp = gsl_interp_alloc(gsl_interp_linear, size);
    gsl_interp_init(interp, tp,r,size);
}
double parameter_T::get(double now){
    double resA;
    gsl_interp_eval_e(interp, tp,r, now, NULL,&resA);
    return resA ;
}

void C2N::events(std::list<cell_base*>::iterator &c_it, std::list<cell_base*> &current_cells, std::list<cell_base*> &old_cells){
    double death = p.d2n();
    double g2nb2n = p.g2nb2n();
    double g2x2nb2n = p.g2x2nb2n();

    double sum = death + g2nb2n + g2x2nb2n ;
    if (sum == 0) 
        return;
    std::uniform_real_distribution<double> dis(0.0, sum);
    double u = dis(*(p.gen));
    if (u<death) 
        ;
    else if (u<death + g2nb2n){
        current_cells.push_front(new C2N(p,p.c_time));   
        current_cells.push_front(new C2N(p,p.c_time));     
    }else{
        current_cells.push_front(new C2x2N(p,p.c_time));
    }
    t_death = p.c_time; 
    old_cells.splice(old_cells.end(),current_cells,c_it);
    return;
    
}


void C2x2N::events(std::list<cell_base*>::iterator &c_it, std::list<cell_base*> &current_cells, std::list<cell_base*> &old_cells){
    double death = p.d2x2n();
    
    double k2nb2x2n = p.k2nb2x2n();
    double k2x2nb2x2n = p.k2x2nb2x2n();
    double k4nb2x2n = p.k4nb2x2n();
    
    double sum = death + k2nb2x2n + k2x2nb2x2n + k4nb2x2n;
    if (sum == 0) 
        return;
    std::uniform_real_distribution<double> dis(0.0, sum);
    double u = dis(*(p.gen));
    if (u<death) 
        ;
    else if (u<death + k2nb2x2n){
        current_cells.push_front(new C2N(p,p.c_time));   
        current_cells.push_front(new C2N(p,p.c_time));     
        current_cells.push_front(new C2N(p,p.c_time));  
        current_cells.push_front(new C2N(p,p.c_time));  
    }else if (u<death + k2nb2x2n + k2x2nb2x2n){
        current_cells.push_front(new C2x2N(p,p.c_time));
        current_cells.push_front(new C2x2N(p,p.c_time));
    }else{
        current_cells.push_front(new C4N(p,p.c_time));
        current_cells.push_front(new C4N(p,p.c_time));
    }
    t_death = p.c_time; 
    old_cells.splice(old_cells.end(),current_cells,c_it);   

    return;
}

void C4N::events(std::list<cell_base*>::iterator &c_it, std::list<cell_base*> &current_cells, std::list<cell_base*> &old_cells){
    double death = p.d4n();
    
    double k2nb4n = p.k2nb4n();
    double k2x2nb4n = p.k2x2nb4n();
    double k4nb4n = p.k4nb4n();
    
    double sum = death + k2nb4n + k2x2nb4n + k4nb4n;
    if (sum == 0) 
        return;
    std::uniform_real_distribution<double> dis(0.0, sum);
    double u = dis(*(p.gen));
    if (u<death) 
        ;
    else if (u<death + k2nb4n){
        current_cells.push_front(new C2N(p,p.c_time));   
        current_cells.push_front(new C2N(p,p.c_time));     
        current_cells.push_front(new C2N(p,p.c_time));  
        current_cells.push_front(new C2N(p,p.c_time));  
    }else if (u<death + k2nb4n + k2x2nb4n){
        current_cells.push_front(new C2x2N(p,p.c_time));
        current_cells.push_front(new C2x2N(p,p.c_time));
    }else{
        current_cells.push_front(new C4N(p,p.c_time));
        current_cells.push_front(new C4N(p,p.c_time));
    }
    t_death = p.c_time; 
    old_cells.splice(old_cells.end(),current_cells,c_it);   

    return;
}


void sim::run( parameters &p, double delta_t,double delta_t_out, double maxT, unsigned int seed, unsigned int num_2n, unsigned int num_2x2n, unsigned int num_4n){
    std::list<cell_base*> current_cells;
    std::list<cell_base*> old_cells;
    
    std::mt19937 gen(seed);
    
    // inizsalizing
    p.set_random(&gen,delta_t);
    p.max_age = maxT;
 

    double OODT = 1/delta_t_out;
    unsigned int t_int = maxT*OODT;
    for (unsigned int j=0;j<NUMCELLS;j++){
        std::vector<std::vector<int>> tmp;
        for (unsigned int i=0;i<t_int;i++){
            tmp.push_back(std::vector<int> (t_int,0));
        }
        time_age.push_back(tmp );
    }
    //init cells
    
    for (unsigned int i=0;i<num_2n;i++){
        current_cells.push_back(new C2N(p,0 ));
    }
    for (unsigned int i=0;i<num_2x2n;i++){
        current_cells.push_back(new C2x2N(p,0 ));
    }
    for (unsigned int i=0;i<num_4n;i++){
        current_cells.push_back(new C4N(p,0 ));
    }
    
    for (double t=0;t<=maxT;t+=delta_t){
        p.get(t);
        auto it = current_cells.begin();
        while (it != current_cells.end()){
            auto it_c = it;
            auto element = *it;
            it++;
            element->events(it_c,current_cells,old_cells);        
        }
        //clear old data and write to output array
        for (auto& oc : old_cells) {
            check_cell(OODT, oc,time_age);
        }
        old_cells.clear();
    }
         
    for (auto& cc : current_cells) {
        check_cell(OODT, cc,time_age);
    }
    old_cells.clear();
  
}

inline void addcell_to_hist(cell_type celltype ,unsigned int time_index_born ,unsigned int time_index_deathmax, std::vector<std::vector<std::vector<int>>>& time_age){
    if (time_index_born == time_index_deathmax){
        time_age[celltype][time_index_born][0] += 0;
    }else{
        for (unsigned int time_index=time_index_born;time_index<time_index_deathmax;time_index++){
            time_age[celltype][time_index][time_index-time_index_born] += 1;
        }
    }
}

void check_cell(double OODT, cell_base *cell, std::vector<std::vector<std::vector<int>>>& time_age){
    unsigned int time_index_born = cell->t_birth*OODT;
    unsigned int time_index_deathmax = cell->t_death*OODT;
    time_index_deathmax = time_index_deathmax < time_age[0].size() ? time_index_deathmax : time_age[0].size();
   
    if (cell->type == t_2n){
        addcell_to_hist(t_2n,time_index_born, time_index_deathmax, time_age);
    }else if  (cell->type == t_2x2n){
        addcell_to_hist(t_2x2n,time_index_born, time_index_deathmax, time_age);
    }else{
        addcell_to_hist(t_4n,time_index_born, time_index_deathmax, time_age);
    }
}