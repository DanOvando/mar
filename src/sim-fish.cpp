#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
NumericVector cpp_seq(int steps,double step_size) {
  
  NumericVector y(steps);
  
  y(0) = 1;
  
  for (int i = 1; i < y.length(); i++){
    
    y(i) = y(i - 1) + step_size;
    
  }
  
  return y ;
}


// [[Rcpp::export]]
List sim_fish(
    const NumericVector length_at_age, 
    const NumericVector weight_at_age,
    const NumericVector maturity_at_age,
    const NumericMatrix f_p_a, // fishing mortality by patch and age
    const List seasonal_movement,
    const List movement_seasons,
    bool dd_movement,
    Rcpp::NumericMatrix last_n_p_a, // last numbers by patch and age
    const int patches,
    const int burn_steps, // number of burn steps if burn is in effect
    const double time_step, // the time step in portions of a year
    double season, // what season you're currently in
    const double steepness,
    const NumericVector r0s,
    double ssb0, // unfished spawning stock biomass across entire system
    NumericVector ssb0_p, // unfished spawning stock biomass in each patch
    const NumericVector m_at_age,
    bool tune_unfished, //0 = use a spawner recruit relationship, 1 = don't
    const int rec_form) // recruitment form, one of ....
  { 

  int ages = length_at_age.length();

  MatrixXd  tmpmat(as<MatrixXd>(last_n_p_a)); // create temporary matrix for movement

  NumericMatrix n_p_a(patches, ages); // numbers in patch p age a

  NumericMatrix b_p_a(patches, ages); // biomass in patch p age a

  NumericMatrix ssb_p_a(patches, ages); // spawning stock biomass in patch p age a

  NumericMatrix c_p_a(patches, ages); // catch at age 

  Eigen::MatrixXd movement = Eigen::MatrixXd::Zero(patches, patches);

  Eigen::MatrixXd density_gradient = Eigen::MatrixXd::Zero(patches, patches);
  
  //////////////////// tune things ////////////////////////
  NumericMatrix tmp_n_p_a = clone(last_n_p_a);

  Rcpp::List tmppop; // annoying step related to memory pointers

  
  // find the correct movement matrix for the season that you're in
  
  for (int s  = 0; s < movement_seasons.length(); s++){ // find which season block you're in
    
    NumericVector tmp_seas = movement_seasons[s];
    
    bool b = Rcpp::any( tmp_seas == season).is_true();
    
    if (b == 1){
      
      movement = seasonal_movement[s];
      
      s = movement_seasons.length() + 1; // stop loop once you've found what season you're in
    } 
    
  } // close seasonal movement finder loop
  
  
  // tune unfished conditions
  
  
  
  if (Rcpp::traits::is_na<REALSXP>(ssb0) == 1){ // basically, if is.na(ssb0), tune unfished conditions

    NumericVector burn_seq = cpp_seq(burn_steps,time_step);
    
    for (int b = 0; b < burn_steps; b++){

      season = burn_seq(b) - floor(burn_seq(b)); // determine what season you're in
      
       
      // these HAVE to be in the same order
      // as function call: get lots of warnings
      // if you try and name them
      // 
      // 
      
      tmppop = sim_fish(
        length_at_age,
        weight_at_age,
        maturity_at_age,
        f_p_a,
        seasonal_movement,
        movement_seasons,
        dd_movement,
        tmp_n_p_a,
        patches,
        0,
        time_step,
        season,
        steepness,
        r0s,
        -999,
        ssb0_p,
        m_at_age,
        1, // this HAS to be 1 inside burn loop
        rec_form);

      tmp_n_p_a = wrap(tmppop["n_p_a"]);

    }

    NumericMatrix tmp_ssb_p_a = tmppop["ssb_p_a"];

    // collect unfished spawning stock biomass 
    ssb0 = sum(tmp_ssb_p_a); 
    
    ssb0_p = rowSums(tmp_ssb_p_a);

  }

  //////////////////// move ////////////////////////


  //// apply density-dependent movement as needed ////
  
  
  // you need a p by p matrix with the relative ssb in each patch relative to the source patch in each column,
  // such that the sum of each column is 1
  
  // compare headroom in each patch
  
  if (tune_unfished == 0 &&  dd_movement == 1){
    
    double ssb_source = 0;
    
    double ssb_dest = 0;
    
    double space_source = 0;
    
    double space_dest = 0;
    
    double tmp_delta;
    
    for (int p = 0; p <patches; p++){
      
      ssb_source = sum(last_n_p_a(p,_) * weight_at_age * maturity_at_age); // calculate ssb in the source patch
      
      tmp_delta = ssb0_p(p) - ssb_source; // calculate headroom in the source patch
      
      space_source = std::max(tmp_delta, 0.0); // make sure headroom is positive, where more headroom less incentive to move
      
      double colsum = 1e-6; // seriously annoying step since I can't figure out how to do colsums of a matrix
      
      double colsum2 = 1e-6;
      
      for (int pp = 0; pp <patches; pp++){
        
        
        ssb_dest= sum(last_n_p_a(pp,_) * weight_at_age * maturity_at_age); // calculate ssb in the destination patch
        
        tmp_delta = ssb0_p(pp) - ssb_dest; // calculate headroom in the destination patch
        
        space_dest = std::max(tmp_delta,0.0);
        
        density_gradient(pp,p) = std::max(space_dest - space_source,0.0); // compare headroom, such that if both places are full, no incentive to move, bigger headroom gradient, more incentive to move
        
        colsum += density_gradient(pp,p);
        
      } // close inner patch loop
      
      // seriously ugly here but works until I can figure out how to subset and sum columns in rcppeigen
      for (int pp = 0; pp<patches; pp++){
        
        density_gradient(pp,p)  =   density_gradient(pp,p) / colsum; // normalize each column to sum to 1
        
        movement(pp,p) = movement(pp,p) +  density_gradient(pp,p); // weighting would go here
          
        colsum2 += movement(pp,p); // calculate colsums of updated movement matrix
          
      }
      
      for (int pp = 0; pp<patches; pp++){
        
        movement(pp,p) = movement(pp,p) /  colsum2; // normalize each column of the updated movement matrix
          
      }
      
    } // close outer patch loop
    

  } // close if dd_movement statement
  
  
  
  // tmpmat =  tmpmat.transpose() * movement; // matrix multiplication of numbers at age by movement matrix

  tmpmat =  movement * tmpmat; // matrix multiplication of numbers at age by movement matrix
  
  SEXP tmp = Rcpp::wrap(tmpmat); // convert from eigen to Rcpp

  // SEXP tmp = Rcpp::wrap(tmpmat.transpose()); // convert from eigen to Rcpp
  
  last_n_p_a = tmp; // set last population to post-movement last population

  //////////////////// grow ////////////////////////


  NumericVector plus_group = last_n_p_a(_,ages - 1) * exp(-(m_at_age(ages - 1) +f_p_a(_,ages - 1))); // calculate numbers in the oldest group that survive
  
  c_p_a(_,ages - 1) =   (f_p_a(_,ages - 1) / (m_at_age(ages - 1) * time_step + f_p_a(_,ages - 1))) * last_n_p_a(_,ages - 1) * (1 - exp(-(m_at_age(ages - 1) * time_step + f_p_a(_,ages - 1))));
  
  // age and die
  for (int a = 1; a < ages; a++){

    n_p_a(_,a) =  last_n_p_a(_,a - 1) * exp(-(m_at_age(a - 1) + f_p_a(_,a - 1)));
    
    c_p_a(_,a - 1) =   (f_p_a(_,a - 1) / (m_at_age(a - 1) * time_step + f_p_a(_,a - 1))) * last_n_p_a(_,a - 1) * (1 - exp(-(m_at_age(a - 1) * time_step + f_p_a(_,a - 1))));

  }

  // add numbers that survived in the oldest group to numbers that grew into the oldest group
  n_p_a(_, ages - 1) = n_p_a(_, ages - 1) + plus_group;

  // calculate biomass and spawning stock biomass at age
  for (int p = 0;p < patches; p++){
    
    c_p_a(p,_) = c_p_a(p,_) * weight_at_age;

    b_p_a(p,_) =  n_p_a(p,_) * weight_at_age;

    ssb_p_a(p,_) =  b_p_a(p,_) * maturity_at_age;

  }

  //////////////////// spawn / recruit ////////////////////////

  double ssb = sum(ssb_p_a); // total spawning stock biomass system wide
  
  NumericVector ssb_p = rowSums(ssb_p_a); // spawning stock biomass by patch
  
  NumericVector recruits(patches); // blank vector for recruits by patch
  
    if (tune_unfished == 1){ // turn off stock recruitment relationship

      if (rec_form == 0){
      
      // recruits = rep(r0 / patches, patches);
      
      recruits = r0s;
      
      } else if (rec_form == 1){
        
        // recruits = rep(r0 / patches, patches);
        
        recruits = r0s;
        
        
      }
      
    } else { // if stock recruitment relationship is in effect


      if (rec_form == 0){  // global beverton-holt density dependence, distribute recruits according to recruitment habitat

      recruits = (((0.8 * sum(r0s) * steepness * ssb) / (0.2 * ssb0 * (1 - steepness) + (steepness - 0.2) * ssb))) * (r0s / sum(r0s));
        
      } else if (rec_form == 1){ // local beverton-holt density dependence, r0 set by local habitat
        
        recruits = ((0.8 * r0s * steepness * ssb_p) / (0.2 * (ssb0_p + 1e-6) * (1 - steepness) + (steepness - 0.2) * (ssb_p + 1e-6)));
        
      } // close recruitment ifs
    }
    
    // assign recruits
    n_p_a(_,0) = recruits;
    
    b_p_a(_,0) =   n_p_a(_,0) * weight_at_age(0);
    
    ssb_p_a(_,0) =   b_p_a(_,0) * maturity_at_age(0);
    

  //////////////////// process results ////////////////////////

  return Rcpp::List::create(
    Rcpp::Named("n_p_a") = n_p_a,
    Rcpp::Named("b_p_a") = b_p_a,
    Rcpp::Named("ssb_p_a") = ssb_p_a,
    Rcpp::Named("ssb0") = ssb0,
    Rcpp::Named("ssb0_p") = ssb0_p,
    Rcpp::Named("tmppop") = tmppop,
    Rcpp::Named("c_p_a") = c_p_a,
    Rcpp::Named("dd_gradient") = density_gradient);
} // close fish model

