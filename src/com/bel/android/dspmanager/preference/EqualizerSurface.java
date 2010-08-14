package com.bel.android.dspmanager.preference;

import android.content.Context;
import android.graphics.Canvas;
import android.graphics.Paint;
import android.graphics.Paint.Style;
import android.util.AttributeSet;
import android.view.SurfaceView;

/**
 * Java support for complex numbers.
 * 
 * @author alankila
 */
class Complex {
	final float re, im;
	
	protected Complex(float re, float im) {
		this.re = re;
		this.im = im;
	}
	
	/**
	 * Length of complex number
	 * 
	 * @return length
	 */
	protected float rho() {
		return (float) Math.sqrt(re * re + im * im);
	}
	
	/**
	 * Argument of complex number
	 * 
	 * @return angle in radians
	 */
	protected float theta() {
		return (float) Math.atan2(im, re);
	}
	
	/**
	 * Complex conjugate
	 * 
	 * @return conjugate
	 */
	protected Complex con() {
		return new Complex(re, -im);
	}

	/**
	 * Complex addition
	 * 
	 * @param other
	 * @return sum
	 */
	protected Complex add(Complex other) {
		return new Complex(re + other.re, im + other.im);
	}
	
	/**
	 * Complex multipply
	 * 
	 * @param other
	 * @return multiplication result
	 */
	protected Complex mul(Complex other) {
		return new Complex(re * other.re - im * other.im, re * other.im + im * other.re);
	}
	
	/**
	 * Complex multiply with real value
	 * 
	 * @param a
	 * @return multiplication result
	 */
	protected Complex mul(float a) {
		return new Complex(re * a, im * a);
	}
	
	/**
	 * Complex division
	 * 
	 * @param other
	 * @return division result
	 */
	protected Complex div(Complex other) {
	    float lengthSquared = other.re * other.re + other.im * other.im;
	    return mul(other.con()).div(lengthSquared);
	}
	
	/**
	 * Complex division with real value
	 * 
	 * @param a
	 * @return division result
	 */
	protected Complex div(float a) {
		return new Complex(re/a, im/a);
	}
}

/**
 * Evaluate transfer functions of biquad filters in direct form 1.
 * 
 * @author alankila
 */
class Biquad {
	float b0, b1, b2, a0, a1, a2;

	protected void setPeakingBandEqualizer(float centerFrequency, float samplingFrequency, float dbGain, float bandwidth) {
        double w0 = 2 * Math.PI * centerFrequency / samplingFrequency;
        double A = Math.pow(10, dbGain/40);
        double alpha = Math.sin(w0)/2 * Math.sinh( Math.log(2)/2 * bandwidth * w0/Math.sin(w0) );

        b0 = (float) (1 + alpha*A);
        b1 = (float) (-2*Math.cos(w0));
        b2 = (float) (1 - alpha*A);
        a0 = (float) (1 + alpha/A);
        a1 = (float) (-2*Math.cos(w0));
        a2 = (float) (1 - alpha/A);
	}
	
	protected void setLowShelf(float centerFrequency, float samplingFrequency, float dbGain, float slope) {
        double w0 = 2 * Math.PI * centerFrequency / samplingFrequency;
        double A = Math.pow(10, dbGain/40);
        double alpha = Math.sin(w0)/2 * Math.sqrt( (A + 1/A)*(1/slope - 1) + 2);

        b0 = (float) (A*((A+1) - (A-1)  *Math.cos(w0) + 2*Math.sqrt(A)*alpha));
        b1 = (float) (2*A*((A-1) - (A+1)*Math.cos(w0)));
        b2 = (float) (A*((A+1) - (A-1)  *Math.cos(w0) - 2*Math.sqrt(A)*alpha));
        a0 = (float) ((A+1) + (A-1)     *Math.cos(w0) + 2*Math.sqrt(A)*alpha);
        a1 = (float) (-2*((A-1) + (A+1) *Math.cos(w0)));
        a2 = (float) ((A+1) + (A-1)     *Math.cos(w0) - 2*Math.sqrt(A)*alpha);
	}
	
	protected void setHighShelf(float centerFrequency, float samplingFrequency, float dbGain, float slope) {
        double w0 = 2 * Math.PI * centerFrequency / samplingFrequency;
        double A = Math.pow(10, dbGain/40);
        double alpha = Math.sin(w0)/2 * Math.sqrt( (A + 1/A)*(1/slope - 1) + 2);

        b0 = (float) (A*((A+1) + (A-1)   *Math.cos(w0) + 2*Math.sqrt(A)*alpha));
        b1 = (float) (-2*A*((A-1) + (A+1)*Math.cos(w0)));
        b2 = (float) (A*((A+1) + (A-1)   *Math.cos(w0) - 2*Math.sqrt(A)*alpha));
        a0 = (float) ((A+1) - (A-1)      *Math.cos(w0) + 2*Math.sqrt(A)*alpha);
        a1 = (float) (2*((A-1) - (A+1)   *Math.cos(w0)));
        a2 = (float) ((A+1) - (A-1)      *Math.cos(w0) - 2*Math.sqrt(A)*alpha);
	}
	
	protected Complex evaluateTransfer(Complex z) {
		Complex nom, den;
		Complex zSquared = z.mul(z);
		
		// b0 + b1 / z + b2 / (z*z)
		nom = new Complex(b0, 0)
			.add(new Complex(b1, 0).div(z))
			.add(new Complex(b2, 0).div(zSquared));

		// a0 + a1 / z + a2 / (z*z)
		den = new Complex(a0, 0)
			.add(new Complex(a1, 0).div(z))
			.add(new Complex(a2, 0).div(zSquared));
		
		return nom.div(den);
	}
}

public class EqualizerSurface extends SurfaceView {
	public static int MIN_FREQ = 20;
	public static int MAX_FREQ = 20000;
	public static int SAMPLING_RATE = 44100;
	public static int MIN_DB = -6;
	public static int MAX_DB = 6;
	
	private final int width;
	private final int height;
	
	float[] levels = new float[5];
	private final Paint white, gray, green, red;
	
	public EqualizerSurface(Context context, AttributeSet attributeSet) {
		super(context, attributeSet);

		/* One day when I'm crazy enough I'll work out how to extract this shit from attributeSet. */
		/* Width can't be taken from the canvas itself passed to onDraw() because android passes
		 * different sized canvases to the method. Yay android. */
		width = 300;
		height = 150;

		/* Also, these fuckers by default disable their own drawing. WTF? */
		setWillNotDraw(false);
		
		white = new Paint();
		white.setColor(0xffffffff);
		white.setStyle(Style.STROKE);

		gray = new Paint();
		gray.setColor(0x88ffffff);
		gray.setStyle(Style.STROKE);

		green = new Paint();
		green.setColor(0x8800ff00);
		green.setStyle(Style.STROKE);

		red = new Paint();
		red.setColor(0x88ff0000);
		red.setStyle(Style.STROKE);
	}

	public void setBand(int i, float value) {
		levels[i] = value;		
		postInvalidate();
	}
	
	@Override
	protected void onDraw(Canvas canvas) {
		/* clear canvas */
		canvas.drawRGB(0, 0, 0);
		
		canvas.drawRect(0, 0, width-1, height-1, white);

		for (int freq = MIN_FREQ; freq < MAX_FREQ;) {
			if (freq < 100) {
				freq += 10;
			} else if (freq < 1000) {
				freq += 100;
			} else if (freq < 10000){ 
				freq += 1000;
			} else {
				freq += 10000;
			}
			
			float x = projectX(freq) * width;
			canvas.drawLine(x, 0, x, height - 1, gray);
		}
		
		for (int db = MIN_DB; db < MAX_DB; db += 3) {
			float y = projectY(db) * height;
			canvas.drawLine(0, y, width - 1, y, gray);
		}
		
		Biquad[] biquads = new Biquad[] {
				new Biquad(),
				new Biquad(),
				new Biquad(),
				new Biquad(),
		};
		
		float gain = (float) Math.pow(10, levels[0] / 20);
		biquads[0].setHighShelf(250f / 2f, SAMPLING_RATE, levels[1] - levels[0], 1f);
		biquads[1].setHighShelf(1000f / 2, SAMPLING_RATE, levels[2] - levels[1], 1f);
		biquads[2].setHighShelf(4000f / 2, SAMPLING_RATE, levels[3] - levels[2], 1f);
		biquads[3].setHighShelf(16000f / 2, SAMPLING_RATE, levels[4] - levels[3], 1f);
		
		/* Now evaluate the tone filter. This is a duplication of the filter design in AudioDSP.
		 * The real filter is integer-based and suffers from approximations; this one is realized
		 * nearly perfectly. */
		float oldx = 0;
		float oldy = 0;
		for (float freq = MIN_FREQ; freq < MAX_FREQ; freq *= 1.2f) {
			float omega = freq / SAMPLING_RATE * (float) Math.PI * 2;
			Complex z = new Complex((float) Math.cos(omega), (float) Math.sin(omega));

			float r1 = z.mul(gain).rho();
			float r2 = biquads[0].evaluateTransfer(z).rho();
			float r3 = biquads[1].evaluateTransfer(z).rho();
			float r4 = biquads[2].evaluateTransfer(z).rho();
			float r5 = biquads[3].evaluateTransfer(z).rho();
			
			float dB = lin2dB(r1 * r2 * r3 * r4 * r5);

			float newx = projectX(freq) * width;
			float newy = projectY(dB) * height;
			
			if (oldx != 0) {
				canvas.drawLine(oldx, oldy, newx, newy, green);
			}
			oldx = newx;
			oldy = newy;
		}
		
		for (int i = 0; i < levels.length; i ++) {
			float freq = 62.5f * (float) Math.pow(4, i);
			float x = projectX(freq) * width;
			float y = projectY(levels[i]) * height;
			canvas.drawLine(x, height/2, x, y, red);
			canvas.drawCircle(x, y, 1, red);
		}
	}

	private float projectX(float freq) {
		double pos = Math.log(freq);
		double minPos = Math.log(MIN_FREQ);
		double maxPos = Math.log(MAX_FREQ);
		pos = (pos - minPos) / (maxPos - minPos);
		return (float) pos;
	}

	private float projectY(float dB) {
		float pos = (dB - MIN_DB) / (MAX_DB - MIN_DB);
		return 1 - pos;
	}

	private float lin2dB(float rho) {
		return rho != 0 ? (float) (Math.log(rho) / Math.log(10) * 20) : -99f;
	}
}
