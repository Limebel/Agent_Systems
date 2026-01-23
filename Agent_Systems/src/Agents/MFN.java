package Agents;
import java.util.*;
import java.io.*;

public class MFN {
    private int m; //number of links in the network
    private int[] W; //number of states every link can be in e.x w[0] = 2 means link can be in states 0, 1 or 2
                    //Those states simulate real life issues e.x broken, half-broken, working correctly
    private double[] C; //Maximum capacity (like maximum flow through point)
    private int[] L; //time needed to traverse a link
    private double[] R; //reliability of link (propability of it working)
    private double[] rho; //fault corelation (simulates faults being related I guess)
    private double[] beta; //Some more complicated propability stuff calculated based on R and rho
    private ArrayList<int[]> MPs = new ArrayList<>(); //Minimal paths (there are many in case of breakdown

    public MFN(Builder builder) {
        this.m = builder.m;
        //W,C,L,R,rho refer to specific links so have to have equal length to number of links
        if (builder.W.length != m || builder.C.length != m ||
                builder.L.length != m || builder.R.length != m ||
                builder.rho.length != m) {
            throw new IllegalArgumentException("Vectors of wrong size.");
        }
        //Propability values are between 0 and 1
        for (double value : builder.R) {
            if (value < 0 || value > 1) {
                throw new IllegalArgumentException("Invalid values in R vector.");
            }
        }
        for (double value : builder.rho) {
            if (value < 0 || value > 1) {
                throw new IllegalArgumentException("Invalid values in R vector.");
            }
        }
        //Calculating beta
        this.beta = new double[m];
        for (int i = 0; i < m; i++) {
            beta[i] = 1 + (builder.rho[i] * (1 - builder.R[i])) / builder.R[i];
        }
        this.W = builder.W;
        this.C = builder.C;
        this.L = builder.L;
        this.R = builder.R;
        this.rho = builder.rho;
        this.MPs = (builder.MPs != null) ? builder.MPs : new ArrayList<>();
    }

    public static class Combinatorial{
        public static int factorial (int n){
            int result = 1;
            for (int i = 1; i <= n; i++) {
                result *= i;
            }
            return result;
        }

        public static int newtonSymbol (int a, int b) {
            if (a >= b){
                return factorial(a)/factorial(b)/factorial(a-b);
            }
            else{
                throw new IllegalArgumentException("Inwalid Newton symbol values.");
            }
        }

        // double factorial
        public static int doubleFactorial (int n){
            int result = 1;
            for (int i = n; i >= 1; i-=2) {
                result *= i;
            }
            return result;
        }
    }

    // Calculates propability of specific state (k) for specific node (i)
    private double probabilityOfStateKCi_Formula1 (int k, int i){
        if (i < 0 || i >= m) throw new IllegalArgumentException("Invalid link index");
        if (k < 0 || k > W[i]) throw new IllegalArgumentException("Invalid capacity state");
        if(k >= 1){
            return 1/beta[i]*Combinatorial.newtonSymbol(W[i],k)*
                    Math.pow(R[i]*beta[i], k)*
                    Math.pow(1 - R[i]*beta[i], W[i]-k);
        }
        else{
            return 1-1/beta[i]*
                    (1-Math.pow(1-R[i]*beta[i],W[i]));
        }
    }

    //How fast we can move some amount of flow(d) through path(P) in SSV(X).
    //This is basic time (formula4) plus how fast we can move demanded capacity through bottleneck
    //Obviously if capacity is 0 we can move any flow and time is infinite
    private double transmissionTimeOfd_Formula3 (int[] P, double d, double[] X){
        double capacityOfP = capacityOfPUnderX_Formula5(P, X);
        if (capacityOfP>0){
            return leadTimeOfP_Formula4(P) + Math.ceil(d / capacityOfP);
        }
        else{
            return Double.POSITIVE_INFINITY;
        }
    }

    // Time it takes to cross a Path(P). It's obviously time of crossing every link on a path
    private int leadTimeOfP_Formula4 (int[] P){
        int sum = 0;
        for (int edge : P) {
            sum += L[edge];
        }
        return sum;
    }

    // Path (P) capacity in SSV(X). Path capacity is the minimum capacity on it
    private double capacityOfPUnderX_Formula5 (int[] P, double[] X){
        double min = X[P[0]];
        for (int edge : P) {
            min = Math.min(min, X[edge]);
        }
        return min;
    }

    // Not used formula, because we do not consider budget
/*    private double transmissionCostOfP_Formula6 (int[] P){
        double totalcost = 0;
        for (int edge : P) {
            totalcost += beta[edge];
        }
        return totalcost;
    }*/

    //Time it takes to move some amount of flow(d) through System with States(SSV)-X.
    //It is minimum time for a specific path
    public double transmissionTimeUnderX_Formula8 (double d, double[] X){
        double minTime = Double.POSITIVE_INFINITY;

        for(int[] P : MPs){
            //double pathCost = transmissionCostOfP_Formula6(P);
            //if (pathCost<=b){
                double time = transmissionTimeOfd_Formula3(P, d, X);
                if(time<minTime){
                    minTime = time;
                }
            //}
        }
        return minTime;
    }

    // void getMPs(String fileName) that reads the file with the file name = filename and creates ArrayList<int[]> MPs
    public void getMPs(String filename) {
        try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] tokens = line.trim().split(",");
                int[] nums = Arrays.stream(tokens)
                        .mapToInt(s->Integer.parseInt(s.trim())-1)
                        .toArray();
                MPs.add(nums);
            }
        } catch (IOException e) {
            System.err.println("Error reading file: " + e.getMessage());
        }
    }

    // and a method just for testing:
    public void printMPs() {
        for (int[] row : MPs) {
            System.out.println(Arrays.toString(row));
        }
    }

    //Propability distribution of every state(k) for every link(m)
    public double[][] arPMFLoop() {
        double[][] arPMF = new double[m][];
        for (int i = 0; i < m; i++) {
            arPMF[i] = new double[W[i] + 1];
            for (int k = 0; k <= W[i]; k++) {
                arPMF[i][k] = probabilityOfStateKCi_Formula1(k, i);
            }
        }
        for (int i = 0; i < m; i++) {
            double sum = 0.0;
            for (int k = 0; k < arPMF[i].length; k++) {
                sum += arPMF[i][k];
            }
        }
        return arPMF;
    }

    //Same distribution as above, but continous
    public double[][] CDF(double[][] arPMF) {
        int m = arPMF.length;
        double[][] arCDF = new double[m][];
        for (int i = 0; i < m; i++) {
            int states = arPMF[i].length;
            arCDF[i] = new double[states];
            double cumulative = 0.0;
            for (int k = 0; k < states; k++) {
                cumulative += arPMF[i][k];
                arCDF[i][k] = cumulative;
            }
        }
        return arCDF;
    }

    // This and normalICDF methods could be used to help calculate worst case NSS, but we have phi, so we don't need them
    // formula: phi(x) = 0.5 + 1/(sqrt(2*pi)) * e^(- x^2 / 2) * [ x + x^3/3 + x^5/(3*5) + ... + x^(2n+1)/((2n+1)!!) + ... ]
    public static double normalCDF(double z) {
        int iMax = 30; // when to stop the infinite loop, can be changed e.g. for more accuracy
        double ans = 0.5;
        double multip = (1 / Math.sqrt(2 * Math.PI));
        multip *= Math.exp(- Math.pow(z, 2) / 2);

        double sum = 0.0;
        for (int i = 1; i < iMax; i+=2) {
            sum += Math.pow(z, i) / Combinatorial.doubleFactorial(i);
        }

        multip *= sum;
        ans += multip;
        return ans;
    }

    static double normaICDF(double u){
        if (u <= 0.0 || u >= 1.0) {
            throw new IllegalArgumentException("u must be in (0,1)");
        }

        double tol = 1e-6;       //precision
        double low = -6.0;
        double high = 6.0;
        double mid = 0.0;

        while ((high - low) > tol) {
            mid = (low + high) / 2.0;
            double cdf = normalCDF(mid);
            if (cdf < u) {
                low = mid;
            } else {
                high = mid;
            }
        }

        return mid;
    }



    //Helper functions for finding worst case NSS
    public static double integrandPhi(double y){
        return Math.exp(- Math.pow(y,2) / 2);
    }
    public static double integralPhi(double z, double minusInfinity){
        double ans = 0.0;
        double step = 0.001;
        for(double x = minusInfinity; x < z; x+=step){
            ans += integrandPhi(x) / Math.sqrt(2*Math.PI) * step;
        }
        return ans;
    }
    public static double phi(double th){
        double minusInfinity = -10.0;
        double z = minusInfinity;
        double step = 0.001;
        while(integralPhi(z, minusInfinity)<th){
            z+= step;
        }
        return z;
    }

    //This function help us decide how many SSVs we need to generate for our results to be valid
    public static int worstCaseNSS(double eps, double delta) {
        double phiVal = phi(1 - delta / 2);
        int nN = (int)Math.ceil(Math.pow(phiVal, 2) / Math.pow(2 * eps, 2));
        return nN;
    }

    // Generating random N number of random SSVs. SSV is System State Vector, so just a random state of a network.
    // e.x. SSV=0,1,2. First link has state 0 etc.
    public double[][] randomSSV(int N, double[][] arCDF) {
        int m = arCDF.length;
        double[][] SSV = new double[N][m];
        Random rand = new Random();

        for (int n = 0; n < N; n++) {
            for (int j = 0; j < m; j++) {
                double u = rand.nextDouble();
                int state = 0;
                // find first index k such that CDF >= u
                for (int k = 0; k < arCDF[j].length; k++) {if (arCDF[j][k] >= u) {
                        state = k;
                        break;
                    }
                }
                SSV[n][j] = state * C[j];
            }
        }
        return SSV;
    }


    public static class Builder{
        private int m;
        private int[] W;
        private double[] C;
        private int[] L;
        private double[] R;
        private double[] rho;
        private double[] beta;
        private ArrayList<int[]> MPs;
        public Builder m(int m) {
            this.m = m;
            return this;
        }
        public Builder W(int[] W) {
            this.W = W;
            return this;
        }
        public Builder C(double[] C) {
            this.C = C;
            return this;
        }
        public Builder L(int[] L) {
            this.L = L;
            return this;
        }
        public Builder R(double[] R) {
            this.R = R;
            return this;
        }
        public Builder rho(double[] rho) {
            this.rho = rho;
            return this;
        }
        public Builder MPs(ArrayList<int[]> MPs) {
            this.MPs = MPs;
            return this;
        }
        public MFN build() {
            return new MFN(this);
        }
    }
    public static Builder builder() {
        return new Builder();
    }
}


