package model;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;

import xmipp.MDLabel;
import xmipp.MetaData;

public class PPData {
	
	private List<Family> families;
	private List<Micrograph> micrographs;
	private static PPData ppdata;
	
	private PPData()
	{
		this.families = new ArrayList<Family>();
		this.micrographs = new ArrayList<Micrograph>();
		loadFamilyData();
		loadMicrographsData();
	}
	
	public static  PPData getInstance()
	{
		if(ppdata == null)
			ppdata = new PPData();
		return ppdata;
	}
	
	public List<Family> getFamilies() {
		return families;
	}
	
	public List<Micrograph> getMicrographs()
	{
		return micrographs;
	}

	
	public void saveFamilyData()
	{
		long id;
		String filename = Family.getOFilename();
		try {
			MetaData md = new MetaData();
			for(Family f: families)
			{
				id = md.addObject();
				md.setValueString(MDLabel.MDL_PICKING_FAMILY, f.getName(), id);
				md.setValueInt(MDLabel.MDL_PICKING_COLOR, f.getColor().getRGB(), id);
				md.setValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, f.getSize(), id);
			}
			md.write("families@" + filename);
		} catch (Exception e) {
			PPConfiguration.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
		
	}
	


	public void loadFamilyData()
	{
		families.clear();
		String filename = Family.getOFilename();
		if(!new File(filename).exists())
		{
			families.add(Family.getDefaultFamily());
			return;
		}
		
		Family family;
		int rgb, size;
		String gname;		
		try {
			MetaData md = new MetaData("families@" + filename);
			long[] ids = md.findObjects();
			for (long id: ids) {				
				gname = md.getValueString(MDLabel.MDL_PICKING_FAMILY, id);
				rgb = md.getValueInt(MDLabel.MDL_PICKING_COLOR, id);
				size = md.getValueInt(MDLabel.MDL_PICKING_PARTICLE_SIZE, id);
				family = new Family(gname, new Color(rgb), size);
				families.add(family);
			}				
			if(families.size() == 0)
				throw new IllegalArgumentException(String.format("No families specified on %s", filename));
		} catch (Exception e) {
			PPConfiguration.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
	}
	
	
	public Family getFamily(String name)
	{
		for (Family f : getFamilies())
			if (f.getName().equalsIgnoreCase(name))
				return f;
		return null;
	}
	
	public boolean existsFamilyName(String name) {
		return getFamily(name)!= null;
	}
	
	
	public void loadMicrographsData()
	{
		String xmd = Micrograph.getIFilename();
		micrographs.clear();
		Micrograph micrograph;
		String ctf, filename;		
		try {
			MetaData md = new MetaData(xmd);
			int count = 1; 
			long[] ids = md.findObjects();
			for (long id: ids) {
				
				filename = PPConfiguration.getMicrographPath(md.getValueString(MDLabel.MDL_IMAGE, id));
				ctf = md.getValueString(MDLabel.MDL_PSD_ENHANCED, id);
				micrograph = new Micrograph(filename, ctf);
				loadParticles(micrograph);
				micrographs.add(micrograph);
				count ++;
			}
			if(micrographs.size() == 0)
				throw new IllegalArgumentException(String.format("No micrographs specified on %s", xmd));
			
		} catch (Exception e) {
			PPConfiguration.getLogger().log(Level.SEVERE, e.getMessage(), e);
		}
		
	}
	
	public void saveParticles(Micrograph micrograph) {
		long id;
		try {
				int count = 0;
				MetaData md;
				String block = null;
				for(Family f: families)
				{
					md = new MetaData();
					for (Particle p: micrograph.getParticles()) {
						if(!f.equals(p.getFamily()))
							continue;
						id = md.addObject();
						md.setValueInt(MDLabel.MDL_XINT, p.getX(), id);
						md.setValueInt(MDLabel.MDL_YINT, p.getY(), id);
						
					}
					block = f.getName() + "@" + micrograph.getOFilename();
					if(count == 0)
						md.write(block);
					else
						md.writeBlock(block);
					count ++;
				}
				
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		
	}

	public void loadParticles(Micrograph micrograph) {
		try {
			int x, y;
			Particle particle;
			MetaData md;
			if(!new File(micrograph.getOFilename()).exists())
				return;
			for(Family f: families)
			{
				md = new MetaData(f.getName() +"@" + micrograph.getOFilename());
				long[] ids = md.findObjects();
				for (long id: ids) {				
					
					x = md.getValueInt(MDLabel.MDL_XINT, id);
					y = md.getValueInt(MDLabel.MDL_YINT, id);
					particle = new Particle(x, y, f, micrograph);
					micrograph.addParticle(particle);
				}
			}
		} catch (Exception e) {
			PPConfiguration.getLogger().log(Level.SEVERE, e.getMessage(), e);
			throw new IllegalArgumentException(e.getMessage());
		}
		
	}

}
