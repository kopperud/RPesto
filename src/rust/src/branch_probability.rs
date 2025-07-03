use itertools::all;

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



#[allow(nonstandard_style)]
impl Gradient for BranchProbabilityMultiState{
    fn gradient(&self, du: &mut Vec<f64>, u: &Vec<f64>, _t: &f64 ) -> (){
        let k = u.len() / 2;
        let r = self.eta / (k as f64 - 1.0);

        let sum_E: f64 = u[0..k].iter().sum();
        let sum_D: f64 = u[(k+1)..(2*k)].iter().sum();

        // E(t)
        for i in 0..k{
            du[i] = self.mu[i] - (self.mu[i] + self.lambda[i] + self.eta) * u[i] + self.lambda[i] * u[i] * u[i] + r * (sum_E - u[i]);
        }

        // D(t)
        for i in 0..k{
            du[k+i] =  - (self.mu[i] + self.lambda[i] + self.eta) * u[k+i] + 2.0 * self.lambda[i] * u[k+i] * u[i] + r * (sum_D - u[k+i]);
        }
    }
}


pub struct ForwardProbability{
    pub lambda: Vec<f64>,
    pub mu: Vec<f64>,
    pub eta: f64,
    pub k: usize,
    pub extinction_probability: MonotonicCubicSpline,
}

impl ForwardProbability{
    pub fn new(
        lambda: Vec<f64>,
        mu: Vec<f64>,
        eta: f64,
        extinction_probability: MonotonicCubicSpline,
    ) -> ForwardProbability{
        let k = lambda.len();

        /*
        let mut e = Vec::new();
        for i in 0..times.len(){
            let mut v = Vec::new();
            for j in 0..k{
                v.push(sol[i][j]);
            }
            e.push(v);
        }


        let mut d = Vec::new();
        for i in 0..times.len(){
            let mut v = Vec::new();
            for j in 0..k{
                v.push(sol[i][k+j]);
            }
            d.push(v);
        }

        println!("asd3");
        */

        let res = ForwardProbability{lambda, mu, eta, k, extinction_probability};
        return res;
    }
}

#[allow(nonstandard_style)]
impl Gradient for ForwardProbability{
    fn gradient(&self, dF: &mut Vec<f64>, F: &Vec<f64>, t: &f64 ) -> (){
        let r = self.eta / (self.k as f64 - 1.0);

        let Et = self.extinction_probability.interpolate(*t);
        let sum_F: f64 = F.iter().sum();

        for i in 0..self.k{
            dF[i] = (self.mu[i] + self.lambda[i] + self.eta) * F[i] - 2.0 * self.lambda[i] * F[i] * Et[i] - r * (sum_F - F[i]);
            //dF[i] = - (self.mu[i] + self.lambda[i] + self.eta) * F[i] + 2.0 * self.lambda[i] * F[i] * Et[i] + r * (sum_F - F[i]);
        }
    }
}

// todo: 
// calculate N(t)
// calculate log Bayes factor



