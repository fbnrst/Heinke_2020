#include <stdio.h>
#include <iostream>
#include <list>
#include "Rm.h"
#include <tuple>
#include <stdexcept>
#include <gsl/gsl_interp.h>
#define SLOW
void parameter::set(std::list<double> timep,std::list<double> rate){
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
double parameter::get(double now){
    double res;
    int e = gsl_interp_eval_e(interp, tp,r, now, NULL,&res);
    return res;
}

void C2N::events(std::list<cell_base*>::iterator &c_it, std::list<cell_base*> &current_cells, std::list<cell_base*> &old_cells){
    #if defined(SLOW)
    bool de = p.r_delta2();
    bool ch = p.r_kappa24();
    bool di = p.r_lambda2();
    
    double sum = de*p.delta2+ch*p.kappa24+di*p.lambda2;
    if (sum == 0) 
        return;
    std::uniform_real_distribution<double> dis(0.0, sum);
    double u = dis(*(p.gen));
    if (u<de*p.delta2) 
        ;
    else if (u<de*p.delta2 + ch*p.kappa24){
        current_cells.push_front(new C4N(p,p.c_time));   
    }else{
        current_cells.push_front(new C2N(p,p.c_time));
        current_cells.push_front(new C2N(p,p.c_time)); 
    }
    t_death = p.c_time; 
    old_cells.splice(old_cells.end(),current_cells,c_it);
    return;
    #else
    if (p.r_delta2()){
        t_death = p.c_time; 
        old_cells.splice(old_cells.end(),current_cells,c_it);
    }else if (p.r_kappa24()){
        t_death = p.c_time; 
        old_cells.splice(old_cells.end(),current_cells,c_it);
        current_cells.push_front(new C4N(p,p.c_time));   
    }else if (p.r_lambda2()){
        t_death = p.c_time; 
        old_cells.splice(old_cells.end(),current_cells,c_it);
        current_cells.push_front(new C2N(p,p.c_time));
        current_cells.push_front(new C2N(p,p.c_time)); 
    }else{}
    #endif
}


void C4N::events(std::list<cell_base*>::iterator &c_it, std::list<cell_base*> &current_cells, std::list<cell_base*> &old_cells){
    #if defined(SLOW)
    bool de = p.r_delta4();
    bool ch =p.r_kappa42();
    bool di = p.r_lambda4();
    
    double sum = de*p.delta4+ch*p.kappa42+di*p.lambda4;
    if (sum == 0) 
        return;
    std::uniform_real_distribution<double> dis(0.0, sum);
    double u = dis(*(p.gen));
    if (u<de*p.delta4) 
        ;
    else if (u<de*p.delta4 + ch*p.kappa42){
        current_cells.push_front(new C2N(p,p.c_time));   
        current_cells.push_front(new C2N(p,p.c_time));   
        current_cells.push_front(new C2N(p,p.c_time));   
        current_cells.push_front(new C2N(p,p.c_time));  
    }else{
        current_cells.push_front(new C4N(p,p.c_time));
        current_cells.push_front(new C4N(p,p.c_time)); 
    }
    t_death = p.c_time; 
    old_cells.splice(old_cells.end(),current_cells,c_it);   

    return;
    #else
    if (p.r_delta4()){
        t_death = p.c_time; 
        old_cells.splice(old_cells.end(),current_cells,c_it);
    }else if (p.r_kappa42()){
        t_death = p.c_time; 
        old_cells.splice(old_cells.end(),current_cells,c_it);
        current_cells.push_front(new C2N(p,p.c_time));   
        current_cells.push_front(new C2N(p,p.c_time));   
    }else if (p.r_lambda4()){
        t_death = p.c_time; 
        old_cells.splice(old_cells.end(),current_cells,c_it);
        current_cells.push_front(new C4N(p,p.c_time));
        current_cells.push_front(new C4N(p,p.c_time)); 
    }else{}
    #endif
}

std::pair<std::vector<std::vector<int>>,std::vector<std::vector<int>>>  test( parameters &p, double delta_t,double delta_t_out, double maxT, unsigned int seed, unsigned int num_2n, unsigned int num_4n){
    std::list<cell_base*> current_cells;
    std::list<cell_base*> old_cells;
    
    //std::random_device rd;
    std::mt19937 gen(seed);
    
    // inizsalizing
    p.set_random(&gen,delta_t);
    p.max_age = maxT;
 
    std::vector<std::vector<int>> n2_time_age;
    std::vector<std::vector<int>> n4_time_age; 
    //#define DT  0.1
    //#define OODT  1/DT
    double DT = delta_t_out;
    double OODT = 1/delta_t_out;
    unsigned int t_int = maxT*OODT+1.5;
    for (unsigned int i=0;i<t_int;i++){
        n2_time_age.push_back(std::vector<int> (t_int,0));
        n4_time_age.push_back(std::vector<int> (t_int,0));
    }
    //init cells
    for (unsigned int i=0;i<num_2n;i++){
        current_cells.push_back(new C2N(p,0));
    }
     for (unsigned int i=0;i<num_4n;i++){
       current_cells.push_back(new C4N(p,0));
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
            check_cell(DT, OODT, oc, n2_time_age, n4_time_age);
            delete oc;
        }
        old_cells.clear();
    }
     
    for (auto& cc : current_cells) {
        check_cell(DT, OODT, cc,n2_time_age,n4_time_age);
        delete cc;
    }
    std::pair<std::vector<std::vector<int>>,std::vector<std::vector<int>>> results(n2_time_age,n4_time_age);
    return results;
}

void check_cell(double DT, double OODT, cell_base *cell, std::vector<std::vector<int>>& n2_time_age,std::vector<std::vector<int>>& n4_time_age){
    unsigned int i_max = cell->t_death*OODT+0.5; // 1.5 for overlapping cell live
    if (i_max > n2_time_age.size())
        i_max = n2_time_age.size();
    if (cell->type == t_2n){
        for (unsigned int i=cell->t_birth*OODT+0.5;i<i_max;i++){
            unsigned int age = (i*DT - cell->t_birth)*OODT+0.5;
            n2_time_age[i][age] += 1; 
        }
    }else{
        for (unsigned int i=cell->t_birth*OODT+0.5;i<i_max;i++){
            unsigned int age = (i*DT - cell->t_birth)*OODT+0.5;
            n4_time_age[i][age] += 1; 
        }
    }
}
