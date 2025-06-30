use itertools::izip;

pub trait Gradient{
    fn gradient(&self, du: &mut Vec<f64>, u: & Vec<f64>, t: &f64);
}


// Two solvers
pub trait Solve{
    fn solve_rk4(    &self, u0: Vec<f64>, t0: f64, t1: f64, dense: bool, n_steps: i32) -> (Vec<f64>, Vec<Vec<f64>>);
    fn solve_dopri45(      &self, u0: Vec<f64>, t0: f64, t1: f64, dense: bool, n_steps_init: i32, error_tolerance: f64) -> (Vec<f64>, Vec<Vec<f64>>);
}

// This is a generic implementation, it implements
// the ODE solver algorithms for all types T that
// have the `Gradient` trait
impl<T> Solve for T
where 
    T: Gradient,
{
    // the Runge-Kutta 4 algorithm
    fn solve_rk4(&self, u0: Vec<f64>, t0: f64, t1: f64, dense: bool, n_steps: i32) -> (Vec<f64>, Vec<Vec<f64>>) {
        let n = u0.len();

        let mut k1 = vec![0.0; n];
        let mut k2 = vec![0.0; n];
        let mut k3 = vec![0.0; n];
        let mut k4 = vec![0.0; n];

        let mut u1 = vec![0.0; n];
        let mut u2 = vec![0.0; n];
        let mut u3 = vec![0.0; n];

        let mut u = u0; 
        let mut t = t0;

        let mut sol: Vec<Vec<f64>> = Vec::new();
        let mut times: Vec<f64> = Vec::new();


        let delta_t = (t1 - t0) / (n_steps as f64);

        for _ in 0..n_steps{
            if dense{
                sol.push(u.clone());
                times.push(t.clone());
            }

            self.gradient(&mut k1, &u, &t);
            for i in 0..n{
                u1[i] = u[i] + k1[i] * delta_t * 0.5;
            }

            self.gradient(&mut k2, &u1, &t);
            for i in 0..n{
                u2[i] = u[i] + k2[i] * delta_t * 0.5;
            }
            //let u2 = u + k2 * delta_t * 0.5;

            self.gradient(&mut k3, &u2, &t);
            for i in 0..n{
                u3[i] = u[i] + k3[i] * delta_t * 0.5;
            }
            //let u3 = u + k3 * delta_t;
            self.gradient(&mut k4, &u3, &t);

            for i in 0..n{
                u[i] += (delta_t / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            }
            //u += (delta_t / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
           
            t += delta_t;
        }

        sol.push(u);
        times.push(t);

        return (times, sol);
    }


    // the Dormand-Prince 45 algorithm with error estimation and adaptive step size
    // Dormand, J.R.; Prince, P.J. (1980). "A family of embedded Runge-Kutta formulae". 
    // Journal of Computational and Applied Mathematics. 6 (1): 19–26. 
    fn solve_dopri45(&self, u0: Vec<f64>, t0: f64, t1: f64, dense: bool, n_steps_init: i32, error_tolerance: f64) -> (Vec<f64>, Vec<Vec<f64>>) {
        let n = u0.len();
        let mut k1 = vec![0.0; n];
        let mut k2 = vec![0.0; n];
        let mut k3 = vec![0.0; n];
        let mut k4 = vec![0.0; n];
        let mut k5 = vec![0.0; n];
        let mut k6 = vec![0.0; n];
        let mut k7 = vec![0.0; n];

        let mut u1 = vec![0.0; n];
        let mut u2 = vec![0.0; n];
        let mut u3 = vec![0.0; n];
        let mut u4 = vec![0.0; n];
        let mut u5 = vec![0.0; n];
        let mut u6 = vec![0.0; n];

        //let tol = 1e-4;
        let mut u = u0; 
        let mut t = t0;
        let mut error_estimate: f64;

        let mut sol: Vec<Vec<f64>> = Vec::new();
        let mut times: Vec<f64> = Vec::new();

        let mut delta_four = vec![0.0; n];
        let mut delta_five = vec![0.0; n];

        let mut delta_t = (t1 - t0) / (n_steps_init as f64);

        let mut number_of_steps = 0;

        let mut c = go(&t, &t1, &delta_t);


        //while t < t1 {
        while c{
            // if overshoot, shorten delta_t
            if overshoot(t, t1, delta_t){
                delta_t = t1 - t;
            }

            //println!("delta_t = {}", delta_t);
            //println!("t = {}", t);
            //println!("t1 = {}", t1);

            self.gradient(&mut k1, &u, &t);

            for i in 0..n{
                u1[i] = u[i] + A21*k1[i] * delta_t;
            }
            self.gradient(&mut k2, &u1, &(t + C2*delta_t));

            for i in 0..n{
                u2[i] = u[i] + (A31*k1[i] + A32*k2[i])* delta_t;
            }
            self.gradient(&mut k3, &u2, &(t + C3*delta_t));

            for i in 0..n{
                u3[i] = u[i] + (A41*k1[i] + A42*k2[i] + A43*k3[i])* delta_t;
            }
            self.gradient(&mut k4, &u3, &(t + C4*delta_t));

            for i in 0..n{
                u4[i] = u[i] + (A51*k1[i] + A52*k2[i] + A53*k3[i] + A54*k4[i])* delta_t;
            }
            self.gradient(&mut k5, &u4, &(t + C5*delta_t));

            for i in 0..n{
                u5[i] = u[i] + (A61*k1[i] + A62*k2[i] + A63*k3[i] + A64*k4[i] + A65*k5[i])* delta_t;
            }
            self.gradient(&mut k6, &u5, &(t + C6*delta_t));

            for i in 0..n{
                u6[i] = u[i] + (A71*k1[i] + A72*k2[i] + A73*k3[i] + A74*k4[i] + A75*k5[i] + A76*k6[i])* delta_t;
            }

            self.gradient(&mut k7, &u6,  &(t + C7*delta_t));

            // fourth and fifth-order solutions
            error_estimate = 0.0;

            for i in 0..n{
                delta_four[i] = delta_t * (B1*k1[i]  + B2*k2[i]  + B3*k3[i]  + B4*k4[i]  + B5*k5[i]  + B6*k6[i] + B7*k7[i]);
                delta_five[i] = delta_t * (A71*k1[i] + A72*k2[i] + A73*k3[i] + A74*k4[i] + A75*k5[i] + A76*k6[i]);
                error_estimate += f64::abs(delta_four[i] - delta_five[i]);
            }
            error_estimate = error_estimate / (n as f64);

            //println!("error_estimate = {}", error_estimate);

            // accept solution
            let mut accept = false;

            let mut any_negatives = false;

            for i in 0..n{
                let current_u = u[i] + delta_five[i]; 
                if current_u < 0.0{
                    any_negatives = true;
                    //break;
                }
            }
            //println!("u = {:?}", u);
            //println!("delta_five = {:?}", delta_five);

            /*
            if any_negatives{
                println!("found negative value at delta_t={}$", delta_t);
            }
            */

            if (error_estimate < error_tolerance) & (!any_negatives) {
                accept = true;
            }

            //println!("any negatives = {}", any_negatives);
            //println!("accept = {}", accept);

            if accept{
                if dense{
                    sol.push(u.clone());
                    times.push(t.clone());
                }

                for i in 0..n{
                    u[i] += delta_five[i]; 
                }
                t += delta_t;

                //println!("accept ODE step");

            // reject, go again with smaller delta t
            }else{
                //println!("delta_t before = {}", delta_t);
                delta_t *= 0.4;
                //println!("error too large, go again");
                //println!("delta_t after shortening = {}", delta_t);
            }
            number_of_steps += 1;

            c = go(&t, &t1, &delta_t);
            //println!("delta_t now = {}", delta_t);
        }
        //println!("number of steps (dopri45): {}", number_of_steps);
        
        if t0 == t1{
            sol.push(u.clone());
            times.push(t.clone());
        }

        sol.push(u);
        times.push(t);



        return (times, sol);
    }
}

fn go(t: &f64, t1: &f64, delta_t: &f64) -> bool{
    let mut res = false;

    //println!("delta_t before: {}", delta_t);

    //if delta_t == &0.0{
    if delta_t > &0.0{
        res = t < t1;
    }else if delta_t < &0.0{
        res = t > t1;
    }

    //println!("delta_t after: {}", delta_t);

    return res;
}

fn overshoot(t: f64, t1: f64, delta_t: f64) -> bool{
    let mut res = false;

    if delta_t > 0.0{
        res = (t + delta_t) > t1;
    }else if delta_t < 0.0{
        res = (t + delta_t) < t1;
    }else{
        panic!("asd");
    }

    return res;
}

// The Butcher tableau for the dopri45 algorithm
const _C1: f64 = 0.0;
const C2: f64 = 1.0/5.0;
const C3: f64 = 3.0/10.0;
const C4: f64 = 4.0/5.0;
const C5: f64 = 8.0/9.0;
const C6: f64 = 1.0;
const C7: f64 = 1.0;

const A21: f64 = 1.0/5.0;

const A31: f64 = 3.0/40.0;
const A32: f64 = 9.0/40.0; 

const A41: f64 = 44.0/45.0;
const A42: f64 = -56.0/15.0;
const A43: f64 = 32.0/9.0;

const A51: f64 = 19372.0/6561.0;
const A52: f64 = -25360.0/2187.0;
const A53: f64 = 64448.0/6561.0;
const A54: f64 = -212.0/729.0;

const A61: f64 = 9017.0/3168.0;
const A62: f64 = -355.0/33.0;
const A63: f64 = 46732.0/5247.0;
const A64: f64 = 49.0/176.0;
const A65: f64 = -5103.0/18656.0;

const A71: f64 = 35.0/384.0;
const A72: f64 = 0.0;
const A73: f64 = 500.0/1113.0;
const A74: f64 = 125.0/192.0;
const A75: f64 = -2187.0/6784.0;
const A76: f64 = 11.0/84.0;

const B1: f64 = 5179.0/57600.0;
const B2: f64 = 0.0;
const B3: f64 = 7571.0/16695.0;
const B4: f64 = 393.0/640.0;
const B5: f64 = -92097.0/339200.0;
const B6: f64 = 187.0/2100.0;
const B7: f64 = 1.0/40.0;


