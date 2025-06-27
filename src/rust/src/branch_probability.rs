use crate::odesolver::*;
use crate::spline::*;
use crate::extinction::*;
use crate::spline::*;

pub struct BranchProbability{
    pub lambda: f64,
    pub mu: f64,
}

impl BranchProbability {
    pub fn new(lambda: f64, mu: f64) -> BranchProbability {

        let res = BranchProbability{lambda, mu};
        return res;
    }
}

impl Gradient for BranchProbability{
    fn gradient(&self, du: &mut Vec<f64>, u: & Vec<f64>, _t: &f64, ) -> (){
        //let et = self.extinction_probability.interpolate(*t)[0];
        
        du[0] = self.mu - (self.mu + self.lambda) * u[0]  + self.lambda * u[0] * u[0];
        du[1] = - (self.mu + self.lambda) * u[1] + 2.0 * self.lambda * u[0] * u[1];
    }
}

pub struct BranchProbabilityMultiState{
    pub lambda: Vec<f64>,
    pub mu: Vec<f64>,
    pub eta: f64,
}

impl BranchProbabilityMultiState{
    pub fn new(lambda: Vec<f64>, mu: Vec<f64>, eta: f64) -> BranchProbabilityMultiState {
        let res = BranchProbabilityMultiState{lambda, mu, eta};
        return res;
    }
}



impl Gradient for BranchProbabilityMultiState{
    fn gradient(&self, du: &mut Vec<f64>, u: &Vec<f64>, _t: &f64 ) -> (){
        let k = u.len() / 2;
        let r = self.eta / (k as f64 - 1.0);

        //let et = self.extinction_probability.interpolate(*t);
        // E(t)
        let sum_E: f64 = u[0..k].iter().sum();


        for i in 0..k{
            du[i] = self.mu[i] - (self.mu[i] + self.lambda[i] + self.eta) * u[i] + self.lambda[i] * u[i] * u[i] + r * (sum_E - u[i]);
        }




        let sum_D: f64 = u[(k+1)..(2*k)].iter().sum();

        for i in 0..k{
            du[k+i] =  - (self.mu[i] + self.lambda[i] + self.eta) * u[k+i] + 2.0 * self.lambda[i] * u[k+i] * u[i] + r * (sum_D - u[k+i]);
        }
    }
}
