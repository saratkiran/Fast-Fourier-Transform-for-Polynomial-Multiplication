/**
 * Created by saratkiran on 3/18/14.
 */
import org.apache.commons.math3.complex.Complex;
import java.util.ArrayList;
import java.util.Arrays;
import java.lang.Math;

public class main {
    public static int length = (int)Math.pow(2,10);
    public static int powers = length/2;

    // A class for complex numbers
    private class Complex
    {
       // the real and imaginary parts
        private double real;
        private double imaginary;

         //Constructs a complex number given the real and imaginary parts.
       public Complex(double real, double imaginary)
        {
            this.real = real;
            this.imaginary = imaginary;
        }

       //  sum of two complex numbers
        public Complex add(Complex c)
        {
            return new Complex(real+c.real, imaginary+c.imaginary);
        }

     // difference of two complex numbers
        public Complex subtract(Complex c)
        {
            return new Complex(real-c.real, imaginary-c.imaginary);
        }

        //product of two complex numbers
        public Complex multiply(Complex c)
        {
            return new Complex(real*c.real - imaginary*c.imaginary,
                    real*c.imaginary + imaginary*c.real);
        }

        public  Complex multiply_3(Complex c){
            double ac = real*c.real;
            double bd = imaginary*c.imaginary;
            return new Complex(ac-bd,((real+imaginary)*(c.real+c.imaginary))-ac-bd);
        }

        // get the real part of complex number
        public double getReal()
        {
            return real;
        }

        // To display the complex number as string
        public String toString(){
            return String.format("%8.4f  + %8.4fi", real, imaginary);
        }
    }

    public static void main(String[] args) {
        for(int l =0;l<1;l++){ // Number of trails
        Complex[] x = new Complex[2*length];

        double[] p = new double[length];  // P polynomial co- efficients
        double[] q = new double[length];    // Q polynomial co- efficients


        // double[] p = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0};  // P polynomial co- efficients
        //double[] q = {1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0};    // Q polynomial co- efficients


        //random double numbers generation
        random(p,-1.0,1.0,length);
        random(q,-1.0,1.0,length);
        main d = new main();
        int n = length*2;
        int logn= (int)(Math.log(n)/Math.log(2)+1e-10);

        int[] shuffle = precomputed_rbs(n,logn);
        Complex[] omega = new Complex[n];
        Complex[] omega_inv = new Complex[n];
        omega = d.omega_computation(omega,n);
        omega_inv = d.omega_inverse_computation(omega_inv,n);
        //System.out.println("omega:");

        long startTime = System.nanoTime();    //Start the timer
         // to get the rescursive solution
        double[] result = d.multiply(p,q,n,omega,omega_inv);
           // to get the dynamic programming solution
        double[] dp_result = d.dp_multiply(p,q,n,omega,shuffle,omega_inv);
        /*System.out.println("three solution: ");
        double [] solution_re3 = three_recursive(p, q, length);
        System.out.println(Arrays.toString(solution_re3));*/

          //displaying the recursive solution
        System.out.println("Normal:");
        System.out.println(Arrays.toString(result));

          // Finding mean absolute error
            /*
        double[] diff = new double[result.length];
        double mea = 0;
        for(int j=0;j<result.length;j++){
            diff[j] = solution_re3[j] - result[j];
            if(diff[j] < 0)
                diff[j] = -1 * diff[j];
            mea += diff[j];
        }
        System.out.println("mean absolute error: ");
        System.out.println(mea/result.length);*/
        //System.out.println(Arrays.toString(diff));

           System.out.println("DP:");
        System.out.println(Arrays.toString(dp_result));

        long stopTime = System.nanoTime();
        float elapsedTime = stopTime - startTime;
        System.out.println(elapsedTime/1000000);    // Calculate the time taken.
    }

    }
    // Function to multiply two numbers.
    public double[] multiply(double[] p, double[] q,int n,Complex[] omega, Complex[] omega_inv)
    {
        //Generate complex objects with 2 times the size
        Complex[] p_double = new Complex[n];
        Complex[] q_double = new Complex[n];

        // Copy the coefficents into the real part of complex number and pad zeros to remainaing terms.
        for (int i=0; i<n/2; i++)
           p_double[i] = new Complex(p[i], 0);
        for (int i=n/2; i<n; i++)
            p_double[i] = new Complex(0, 0);
        for (int i=0; i<n/2; i++)
            q_double[i] = new Complex(q[i], 0);
        for (int i=n/2; i<n; i++)
            q_double[i] = new Complex(0, 0);

        int pow = 1;

        // Apply the FFT to the two factors
        Complex[] solp = fft(p_double, omega, n, pow);
        Complex[] solq = fft(q_double, omega, n,pow);

        // Multiply the results pointwise recursive
        Complex[] finalsol = new Complex[n];
        for (int i=0; i<n; i++)
            finalsol[i] = solp[i].multiply(solq[i]);

        // Apply the FFT to the pointwise product
        Complex[] poly = fft(finalsol, omega_inv, n,pow);

        //
        //   get the results by normalising the results
        //
        double[] result = new double[n-1];
        for (int i=0; i<n-1; i++)
            result[i] = poly[i].getReal()/n;

        /*result[0] = poly[0].getReal()/n;
        for (int i=1; i<n-1; i++)
            result[i] = poly[n-i].getReal()/n;*/

        //
        //   get the results by normalising the results.
        //

        return result;
    }

    public double[] dp_multiply(double[] p, double[] q,int n,Complex[] omega,int[] shuffle,Complex[] omega_inv)
    {
        //Generate complex objects with 2 times the size
        Complex[] p_double = new Complex[n];
        Complex[] q_double = new Complex[n];

        // Copy the coefficents into the real part of complex number and pad zeros to remainaing terms.
        for (int i=0; i<n/2; i++)
            p_double[i] = new Complex(p[i], 0);
        for (int i=n/2; i<n; i++)
            p_double[i] = new Complex(0, 0);
        for (int i=0; i<n/2; i++)
            q_double[i] = new Complex(q[i], 0);
        for (int i=n/2; i<n; i++)
            q_double[i] = new Complex(0, 0);

        // Apply the Dynamic Programming FFT
        Complex[] dp_solp = dpfft(p_double, omega, n,shuffle);
        Complex[] dp_solq = dpfft(q_double, omega, n,shuffle);

        // Multiply the results pointwise for dynamic Programming fft
        Complex[] dp_finalsol = new Complex[n];
        for (int i=0; i<n; i++)
            dp_finalsol[i] = dp_solp[i].multiply_3(dp_solq[i]);

        Complex[] dp_poly = dpfft(dp_finalsol, omega_inv, n, shuffle);

        //
        //   get the results by normalising the results for dynamic programming
        //


        double[] dp_result = new double[n-1];
        for (int i=0; i<n-1; i++)
            dp_result[i] = dp_poly[i].getReal()/n;
       /* dp_result[0] = dp_poly[0].getReal()/n;
        for (int i=1; i<n-1; i++)
            dp_result[i] = dp_poly[n-i].getReal()/n;
*/
        return dp_result;
    }

    // Function to obtain fft
    private Complex[] fft(Complex[] pol, Complex[] omega, int n,int pow){

        // the base case --
        if (n==1)
            return pol;

        // split  coefficients into two pieces - even and idd
        Complex[] poly_even = new Complex[n];
        Complex[] poly_odd = new Complex[n];

        // send them into two different array of objects
        for(int k =0;k<n;k++){
            if(k%2 ==0)
                poly_even[k/2] = pol[k];
            else
                poly_odd[k/2] = pol[k];
        }

        // get the squares of omega values
        /*Complex[] xsquare = new Complex[n/2];
        for(int i=0;i<n/2;i++)
           xsquare[i]=omega[i].multiply(omega[i]);*/

        // call the even and odd objects arrays recursively
        Complex[] solution_even = fft(poly_even, omega, n/2,pow*2);
        Complex[] solution_odd =  fft(poly_odd, omega, n/2,pow*2);

        // create space for the results
        Complex[] poly_sol = new Complex[n];
        for (int i=0; i<n; i++)
           poly_sol[i] = new Complex(0,0);

        // combine the pieces
        for (int i=0; i<n/2; i++)
        {
            poly_sol[i] = solution_even[i].add(omega[i*pow].multiply(solution_odd[i]));
            poly_sol[i+n/2] = solution_even[i].subtract(omega[i*pow].multiply(solution_odd[i]));

        }

        // return the final solution of fft
        return poly_sol;
    }


    // Pre compute the omega values
    public Complex[] omega_computation(Complex[] x,int length){
        for(int i=0;i<length;i++){
            x[i] = new Complex(Math.cos((2*i*Math.PI)/length),Math.sin((2*i*Math.PI)/length));
        }
        return x;
    }

    // Pre compute the omega values
    public Complex[] omega_inverse_computation(Complex[] x,int length){
        for(int i=0;i<length;i++){
            x[i] = new Complex(Math.cos((2*i*Math.PI)/length),-1*Math.sin((2*i*Math.PI)/length));
        }
        return x;
    }



    // Random values generation
    public static void random(double[] p,double min, double max,int length ){
        double diff = max - min;
        for(int i=0;i<length;i++){
            p[i] = min + Math.random( ) * diff;
        }
    }

    // NOT USING THIS FUNCTION
    public static Complex[] recfft(Complex[] p_double,Complex[] x,int len){
        //base case

        if(len == 1){
            Complex[] poly = new Complex[1];
                  poly[0]  = p_double[0];
            return poly;
        }
        else{
            Complex[] poly = new Complex[len];
            Complex[] p_even = new Complex[len/2];
            Complex[] p_odd = new Complex[len/2];
            Complex[] xsquare = new Complex[len/2];
            for(int k =0;k<len;k++){
                if(k%2 ==0)
                    p_even[k/2] = p_double[k];
                else
                    p_odd[k/2] = p_double[k];
            }
            System.out.println(Arrays.toString(p_even));
            System.out.println(Arrays.toString(p_odd));
            for(int k =0;k <len/2;k++){
                xsquare[k]= x[k].multiply(x[k]);
            }
            Complex[] poly_e = recfft(p_even, xsquare, len / 2);
            System.out.println(Arrays.toString(poly_e));
            Complex[] poly_o = recfft(p_odd, xsquare, len / 2);
            System.out.println(Arrays.toString(poly_o));

            for(int k= 0;k<len/2;k++){

                poly[k] = poly_e[k].add(poly_o[k].multiply(x[k*powers]));
                poly[k+(len/2)] = poly_e[k].subtract(poly_o[k].multiply(x[k*powers]));

            }
            powers = powers/2;
            return poly;
        }
    }

    // Dynamic FFT
    public Complex[] dpfft(Complex[] poly, Complex[] x,int n,int[] shuffle){
        int logn= (int)(Math.log(n)/Math.log(2)+1e-10);

        //Complex[][] sol = new Complex[logn+1][n];
        Complex[][] sol = new Complex[2][n];
        for(int i=0; i<n; i++){
            sol[0][shuffle[i]] = poly[i];
        }
        int power = n/2;
        int size = 2;
        int k;

        for(k=1;k<=logn;k++){
            for(int i =0;i<n ;i=i+size){
                for(int j=0;j<size/2;j++){
                    Complex odd = x[j*power].multiply_3(sol[(k-1)%2][(i+j+size/2)]);
                    sol[k%2][i+j] = sol[(k-1)%2][i+j].add(odd);
                    sol[k%2][i+j+size/2] = sol[(k-1)%2][i+j].subtract(odd);
                }
            }
            power = power/2;
            size = size*2;
        }
        Complex[] temp = new Complex[n];
        for(int i =0;i <n;i++)
            temp[i] = sol[(k-1)%2][i];

        return temp;

    }

    // Calculate  the precomputed RBS value.
    public static int[] precomputed_rbs(int length,int logn){

        int[] shuffled = new int[length];
        for(int i=0; i< length;i++)
            shuffled[i] = RBS(i,logn);
        return shuffled;
    }

    // Calculate the RBS values
    public static int RBS(int i,int k){
        if(k == 0 )
            return i;
        if (i%2 ==1)
            return (int)Math.pow(2,k-1)+ RBS(i/2,k-1);
        else
            return RBS(i/2,k-1);
    }

    //LAST ASSIGMNMWNENT
    public static double[] three_recursive(double[] p,double[] q,int length){

        //base case
        if(length == 1){
            double[] pq = new double[1];
            pq[0] = p[0]*q[0];
            return pq;
        }
        else{
            //array declarations for pl,ph,ql,qh
            double[] pq = new double[(2*length)];
            double[] pl = Arrays.copyOfRange(p, 0, (length/2));
            double[] ph = Arrays.copyOfRange(p, (length/2), length);
            //System.out.println(Arrays.toString(ph));
            double[] ql = Arrays.copyOfRange(q, 0, (length/2));
            double[] qh = Arrays.copyOfRange(q, (length/2), length);

            //recursive calls
            double[] plql = three_recursive(pl, ql, length / 2);
            double[] phqh = three_recursive(ph, qh, length / 2);

            // pl+ph , ql+qh
            double[] plandph = new double[length];
            double[] qlandqh = new double[length];

            for(int i=0;i<length/2;i++){
                plandph[i] = pl[i] + ph[i];
                qlandqh[i] = ql[i] + qh[i];
            }

            //recursive call - (pl+ph)(ql+qh)
            double[] plmulqh = three_recursive(plandph,qlandqh,length/2);

            // (pl+qh)(ql+qh) - plql - phqh
            double[] plmulqh_f = new double[(length)];
            for (int k =0;k<length-1; k++){
                plmulqh_f[k] = plmulqh[k] -plql[k] -phqh[k];
            }

            // Adding the solutions as plql + ((pl+ph)(ql+qh)-plql-qlqh)x^n/2 +phqhx^n

            int mid = ((length));
            int x=1;
            while(x<=((length*2)-1)){
                if(x <= length/2) {
                    pq[x] = plql[x-1];
                }
                else if(x >= ((length/2)+1)&& x < length) {
                    pq[x] = plql[x-1] + plmulqh_f[x - ((length / 2)+1)];
                }
                else if (x == mid){
                    pq[x] = plmulqh_f[x-((length/2)+1)];
                }

                else if (x > mid && x <= ((3*(length)/2)-1)) {
                    pq[x] = phqh[x - (mid+1)] + plmulqh_f[x - ((length/2)+1)];
                }
                else {
                    pq[x] = phqh[x - (mid+1)];
                }
                x++;
            }
            // System.out.println(Arrays.toString(pq));
            return Arrays.copyOfRange(pq, 1, pq.length);

        }


    }
}
