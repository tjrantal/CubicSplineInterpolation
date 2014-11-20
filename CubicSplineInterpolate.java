/*
	A class for cubic spline interpolation. Uses JAMA to solve the group of equations
	Written by Timo Rantalainen 2011-2014 tjrantal at gmail dot com . Modified from the one originally published in density distribution analysis source code at densitydistribution.comli.com. Cannot remember which sources I used for the maths, but likely something like http://en.wikipedia.org/wiki/Spline_interpolation and http://mathworld.wolfram.com/CubicSpline.html
	
	GPL 3.0 or above license applies (no copy of the license provided, look at http://www.gnu.org/copyleft/gpl.html for the license terms)
	
	Depends on the Jama package (http://math.nist.gov/javanumerics/jama/). Copy it to a subfolder called Jama
	Create a subfolder build, and compile
	javac -cp "." -d build CubicSplineInterpolate.java
	
	execute (replace .:build with .;build in windows)
	java -cp ".:build" CubicSplineInterpolate
	
	javac -cp "." -d build CubicSplineInterpolate.java && java -cp ".:build" CubicSplineInterpolate
	
*/
import Jama.*;
public class CubicSplineInterpolate{

	/*Main for testing*/
	public static void main(String[] a){
		double[] samplingInstants = new double[11];
		double[] sampledValues = new double[11];
		double[] interpolatedSamplingInstants = new double[10];
		double[] knownValues = new double[10];
		double sigFreq = 15;
		for (int i = 1; i< samplingInstants.length;++i){
			samplingInstants[i] = 0.01*((double)i)+Math.random()*0.01-0.0049;
			sampledValues[i] = Math.sin(2d*Math.PI*sigFreq*samplingInstants[i]);
			if (i<interpolatedSamplingInstants.length){
				interpolatedSamplingInstants[i] = 0.01*((double)i);
				knownValues[i] = Math.sin(2d*Math.PI*sigFreq*interpolatedSamplingInstants[i]);
			}
		}
		double[] interpolatedValues = cubicSplineInterpolate(samplingInstants,sampledValues,interpolatedSamplingInstants);
		/*Print out the result*/
		for (int i = 0; i< interpolatedValues.length;++i){
			System.out.println("si\t"+String.format("%.3f",samplingInstants[i])+"\tKnown\t"+String.format("%.2f",sampledValues[i])+"\tsi\t"+String.format("%.3f",interpolatedSamplingInstants[i])+"\tKnown\t"+String.format("%.2f",knownValues[i])+"\tinterp\t"+String.format("%.2f",interpolatedValues[i]));
		}
		
	}	
	
	/**
		Method for cubic spline interpolation. Cannot remember which sources I used for the maths, but likely likely something like  http://en.wikipedia.org/wiki/Spline_interpolation and http://mathworld.wolfram.com/CubicSpline.html
		@param samplingInstants, a 1D array of sampling instants
		@param sampledValues, a 1D array of sampled values corresponding to the sampling instants
		@param interpolationInstants, a 1D array of time instants to interpolate the sample at. Must be entirely within the samplingInstants.
		@return interpolatedSamples, the samples interpolated at interpolationInstants
	*/
	public static double[] cubicSplineInterpolate(double[] samplingInstants,double[] sampledValues,  double[] interpolationInstants){
		double[][] a = new double[sampledValues.length][sampledValues.length];
		double[] b = new double[sampledValues.length];
		double[] bb = new double[sampledValues.length-1];
		double [] d = new double[sampledValues.length-1];

		a[0][0] = 1;
		a[sampledValues.length-1][sampledValues.length-1] = 1;
		b[0] = 0;
		b[sampledValues.length-1] = 0;
		for (int i = 1; i< sampledValues.length-1;++i){
			a[i][i-1] =samplingInstants[i]-samplingInstants[i-1];
			a[i][i]   =2.0*((samplingInstants[i]-samplingInstants[i-1])+(samplingInstants[i+1]-samplingInstants[i]));
			a[i][i+1] =samplingInstants[i+1]-samplingInstants[i];
			b[i] = (3/(samplingInstants[i+1]-samplingInstants[i]))*(sampledValues[i+1]-sampledValues[i])-(3/(samplingInstants[i]-samplingInstants[i-1]))*(sampledValues[i]- sampledValues[i-1]);
		}
		Matrix A = new Matrix(a);
		Matrix B = new Matrix(b,sampledValues.length);
		Matrix c = A.solve(B);
		for (int i = 1; i < sampledValues.length;i++){
			bb[i-1] = 1/(samplingInstants[i]-samplingInstants[i-1])*(sampledValues[i]-sampledValues[i-1])-(samplingInstants[i]-samplingInstants[i-1])/3*(2*c.get(i-1,0)+c.get(i,0));
			d[i-1] = (c.get(i,0)-c.get(i-1,0))/(3*(samplingInstants[i]-samplingInstants[i-1]));
		}
		double[] interpolatedSamples = new double[interpolationInstants.length];
		double inter;
		/*Reconstruc the interpolated signal*/
		int j = 0;
		for (int i = 0;i<interpolationInstants.length;++i){
			while (samplingInstants[j+1]<interpolationInstants[i]){
				++j;
			}
			inter = (interpolationInstants[i]-samplingInstants[j]);
			interpolatedSamples[i] = sampledValues[j]+bb[j]*inter+c.get(j,0)*Math.pow(inter,2.0)+d[j]*Math.pow(inter,3.0);
		}
		return interpolatedSamples;
	}
}
