use crate::odesolver::*;

pub struct Extinction{
    pub lambda: f64,
    pub mu: f64,
}

impl Gradient for Extinction{
    fn gradient(&self, du: &mut Vec<f64>, u: & Vec<f64>, _t: &f64, ) -> (){
        du[0] = self.mu - (self.mu + self.lambda) * u[0]  + self.lambda * u[0] * u[0];
    }
}

pub struct ExtinctionMultiState{
    pub lambda: Vec<f64>,
    pub mu: Vec<f64>,
    pub eta: f64,
    pub extinction_approximation: bool,
}

impl Gradient for ExtinctionMultiState{
    fn gradient(&self, du: &mut Vec<f64>, u: &Vec<f64>, _t: &f64 ) -> (){
        let k = u.len();
        let r = self.eta / (k as f64 - 1.0);

        let sum_u: f64 = u.iter().sum();

        if self.extinction_approximation{
            for i in 0..k{
                du[i] = self.mu[i] - (self.mu[i] + self.lambda[i]) * u[i] + self.lambda[i] * u[i] * u[i] 
            }
        }else{
            for i in 0..k{
                du[i] = self.mu[i] - (self.mu[i] + self.lambda[i] + self.eta) * u[i] + self.lambda[i] * u[i] * u[i] + r * (sum_u - u[i]);
            }
        }


        /*
        for i in 0..k{
            du[i] = self.mu[i] - (self.mu[i] + self.lambda[i] + self.eta) * u[i] + self.lambda[i] * u[i] * u[i];
        }

        for i in 0..k{
            for j in 0..k{
                du[i] += r * u[j];
            }
        }

        for i in 0..k{
            du[i] -= r * u[i];
        }
        */
    }
}
