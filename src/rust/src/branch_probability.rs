use crate::odesolver::*;
use crate::spline::*;
use crate::extinction::*;
use crate::spline::*;

pub struct BranchProbability{
    pub lambda: f64,
    pub mu: f64,
    pub extinction_probability: MonotonicCubicSpline,
}

impl BranchProbability {
    pub fn new(lambda: f64, mu: f64, height: f64, tol: f64) -> BranchProbability {

        let ode = Extinction{lambda, mu};

        let u0 = vec![0.0];
        let t0 = 0.0;
        let t1 = height;
        let dense = true;
        let n_steps_init = 10;

        let (times, sol) = ode.solve_dopri45(u0, t0, t1, dense, n_steps_init, tol);

        let extinction_probability = MonotonicCubicSpline::new(times, sol, 1);

        let res = BranchProbability{lambda, mu, extinction_probability};

        return res;
    }
}

impl Gradient for BranchProbability{
    fn gradient(&self, du: &mut Vec<f64>, u: & Vec<f64>, t: &f64, ) -> (){
        let et = self.extinction_probability.interpolate(*t)[0];

        du[0] = - (self.mu + self.lambda) * u[0] + 2.0 * self.lambda * u[0] * et;
    }
}

pub struct BranchProbabilityMultiState{
    pub lambda: Vec<f64>,
    pub mu: Vec<f64>,
    pub eta: f64,
}

impl Gradient for BranchProbabilityMultiState{
    fn gradient(&self, du: &mut Vec<f64>, u: &Vec<f64>, _t: &f64 ) -> (){
        let k = u.len();
        let r = self.eta / (k as f64 - 1.0);

        let sum_u: f64 = u.iter().sum();

        for i in 0..k{
            du[i] = self.mu[i] - (self.mu[i] + self.lambda[i] + self.eta) * u[i] + self.lambda[i] * u[i] * u[i] + r * (sum_u - u[i]);
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
