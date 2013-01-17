package xmipp.viewer.particlepicker.extract;

import java.awt.Color;

import xmipp.jni.MDLabel;
import xmipp.viewer.particlepicker.ColorBy;
import xmipp.viewer.particlepicker.Micrograph;
import xmipp.viewer.particlepicker.PickerParticle;

public class ExtractParticle extends PickerParticle
{
	private boolean enabled;
	private long id;
	private double zscore, zscore_shape, zscore_snr1, zscore_snr2, zscore_hist;

	public ExtractParticle(long id, int x, int y, Micrograph m, boolean enabled, double zscore, double zscore_shape, double zscore_snr1, double zscore_snr2, double zscore_hist)
	{
		super(x, y, m);
		this.enabled = enabled;
		this.id = id;
		this.zscore = zscore;
		this.zscore_hist = zscore_hist;
		this.zscore_shape = zscore_shape;
		this.zscore_snr1 = zscore_snr1;
		this.zscore_snr2 = zscore_snr2;
		
	}
	
	public long getId()
	{
		return id;
	}

	public ExtractParticle(long id, int x, int y, Micrograph m, double zscore, double zscore_shape, double zscore_snr1, double zscore_snr2, double zscore_hist)
	{
		this(id, x, y, m, true, zscore, zscore_shape, zscore_snr1, zscore_snr2, zscore_hist);
	}

	public boolean isEnabled()
	{
		return enabled;
	}

	

	public double getZscore()
	{
		return zscore;
	}

	

	public double getZscore_shape()
	{
		return zscore_shape;
	}

	

	public double getZscore_snr1()
	{
		return zscore_snr1;
	}

	
	public double getZscore_snr2()
	{
		return zscore_snr2;
	}

	

	public double getZscore_hist()
	{
		return zscore_hist;
	}

	public void setEnabled(boolean b)
	{
		enabled = b;
		
	}

	public Color getColor(ColorBy colorBy)
	{
		return getColor(colorBy.getId(), colorBy.getColor());
		
	}
	
	

	private double getScore(int id)
	{
		if(id == MDLabel.MDL_ZSCORE)
			return getZscore();
		if(id == MDLabel.MDL_ZSCORE_SHAPE)
			return getZscore_shape();
		if(id == MDLabel.MDL_ZSCORE_SNR1)
			return getZscore_snr1();
		if(id == MDLabel.MDL_ZSCORE_SNR2)
			return getZscore_snr2();
		if(id == MDLabel.MDL_ZSCORE_HISTOGRAM)
			return getZscore_hist();
		return 0;
	}

	public Color getColor(int id, Color color)
	{
		double score = getScore(id)/2 + 0.1 ; //scaled;
		return ColorBy.getColor(score, color);
	}

	
	
	
	

}
