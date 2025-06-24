/* 
 
 modified from 
 https://gitlab.com/ruivieira/mentat/

 */


pub struct MonotonicCubicSpline {
    // vector of knots (time points)
    m_times: Vec<f64>,
    // vector of explicit solutions to the ODE 
    m_values: Vec<Vec<f64>>,
    // something to do with hermite basis function
    m_m: Vec<Vec<f64>>,
    categories: usize
}

impl MonotonicCubicSpline {

    pub fn new(times : Vec<f64>, values : Vec<Vec<f64>>, k: usize) -> MonotonicCubicSpline {

        assert!(times.len() == values.len() && times.len() >= 2 && values.len() >= 2, "Must have at least 2 control points.");

        let n = times.len();

        //let mut secants = vec![0.0 ; n - 1];
        let mut secants = vec![vec![0.0 ; k] ; n - 1];

        //let mut slopes  = vec![0.0 ; n];
        let mut slopes = vec![vec![0.0 ; k] ; n];

        for j in 0..k{
            for i in 0..(n-1) {
                let h = *times.get(i + 1).unwrap() - *times.get(i).unwrap();
                assert!(h > 0.0, "Control points must be monotonically increasing.");
                let vt = values.get(i + 1).unwrap();
                let v = *vt.get(j).unwrap();

                let vti = values.get(i).unwrap();
                let vi = *vti.get(j).unwrap();
                secants[i][j] = (v - vi) / h;
                //*values.get(i + 1).get(j).unwrap() - *values.get(i).unwrap()) / h;
            }
        }
        

        for j in 0..k{
            slopes[0][j] = secants[0][j];

            for i in 1..(n-1) {
                slopes[i][j] = (secants[i - 1][j] + secants[i][j]) * 0.5;
            }
            slopes[n - 1][j] = secants[n - 2][j];
        }

        for j in 0..k{
            for i in 0..(n-1) {
                if secants[i][j] == 0.0 {
                    slopes[i][j] = 0.0;
                    slopes[i + 1][j] = 0.0;
                } else {
                    let alpha = slopes[i][j] / secants[i][j];
                    let beta = slopes[i + 1][j] / secants[i][j];
                    let h = alpha.hypot(beta);
                    if h > 9.0 {
                        let z = 3.0 / h;
                        slopes[i][j] = z * alpha * secants[i][j];
                        slopes[i + 1][j] = z * beta * secants[i][j];
                    }
                }
            }
        }


        let spline = MonotonicCubicSpline {
            m_times: times.clone(),
            m_values: values.clone(),
            m_m: slopes,
            categories: k,
        };
        return spline;
    }

    pub fn hermite(point: f64, time : (f64, f64), value: (&Vec<f64>, &Vec<f64>), m: (&Vec<f64>, &Vec<f64>), k: usize) -> Vec<f64> {
        let h: f64 = time.1 - time.0;
        let z = (point - time.0) / h;

        let mut res = Vec::new();

        for j in 0..k{
            let v = (value.0[j] * (1.0 + 2.0 * z) + h * m.0[j] * z) * (1.0 - z) * (1.0 - z);
            res.push(v);
        }

        return res;
    }

    pub fn interpolate(&self, time : f64) -> Vec<f64> {
        let n = self.m_times.len();

        let p = *self.m_times.get(0).unwrap();
        if time <= p {
            let r = self.m_values.get(0).unwrap();
            return r.clone();
        }
        if time >= *self.m_times.get(n - 1).unwrap() {
            let r = self.m_values.get(n - 1).unwrap();
            return r.clone();
        }

        let mut i = 0;
        while time >= *self.m_times.get(i + 1).unwrap() {
            i += 1;
            if time == *self.m_times.get(i).unwrap() {
                let r = self.m_values.get(i).unwrap();
                return r.clone();
            }
        }

        let t = *self.m_times.get(i).unwrap();
        let tp = *self.m_times.get(i+1).unwrap();
        let timespan = (t, tp);

        let val = self.m_values.get(i).unwrap();
        let valp = self.m_values.get(i+1).unwrap();
        let value = (val, valp);

        let mm = self.m_m.get(i).unwrap();
        let mp = self.m_m.get(i+1).unwrap();
        let m = (mm, mp);

        return MonotonicCubicSpline::hermite(time,
                                                timespan,
                                             value,
                                             m,
                                             self.categories);
    }

}

