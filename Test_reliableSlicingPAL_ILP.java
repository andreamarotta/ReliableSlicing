package com.net2plan.examples.general.offline.nfv;

import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
//import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tint.IntMatrix1D;
import cern.colt.matrix.tint.IntMatrix2D;
import cern.colt.matrix.tint.IntFactory1D;
import cern.colt.matrix.tint.IntFactory2D;

import com.jom.DoubleMatrixND;
import com.jom.Expression;
import com.jom.OptimizationProblem;
import com.net2plan.examples.general.offline.Offline_ipOverWdm_routingSpectrumAndModulationAssignmentILPNotGrooming;
import com.net2plan.interfaces.networkDesign.*;
import com.net2plan.utils.Constants.RoutingType;
import com.net2plan.utils.InputParameter;
import com.net2plan.utils.Pair;
import com.net2plan.utils.Triple;
import org.apache.commons.lang3.StringUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.stream.Collectors;
import java.sql.Timestamp;
import java.text.SimpleDateFormat;

// ASSUMES YOU CAN NOT SPLIT ONE TRAFFIC REQUEST OVER DIFFERENT WAVELENGHTS
// DOES IT MAKES SENSE?

/**
 *
 * 
 * @net2plan.keywords JOM, NFV
 * @net2plan.inputParameters
 * @author Andrea Marotta
 */
public class Test_reliableSlicingPAL_ILP  {



	public static void main(String[] args) throws IOException {
		Locale.setDefault(Locale.US);
		//NetPlan netPlan = NetPlan.loadFromFile(new File("resil.n2p"));
		//NetPlan netPlan = NetPlan.loadFromFile(new File("eon_N18_E66_withTraffic.n2p"));
		//		NetPlan netPlan = NetPlan.loadFromFile(new File("N7.n2p"));/

		NetPlan netPlan = NetPlan.loadFromFile(new File("N7.n2p"));
		Map<String, String> parameters = new HashMap<>();
		parameters.put("solverName", "cplex");
		parameters.put("solverLibraryName","C:\\Program Files\\IBM\\ILOG\\CPLEX_Studio129\\opl\\bin\\x64_win64\\cplex1290.dll");
		parameters.put("maxSolverTimeInSeconds", String.valueOf((double) -1));
		//parameters.put("solverName", "glpk" );
		//parameters.put("solverLibraryName","C:\\WINDOWS\\system32\\glpk_4_48.dll");

		parameters.put("C", String.valueOf((double) 40));
		parameters.put("W", String.valueOf((int) 80));
		//parameters.put("isolation", String.valueOf((int) 4));

		//parameters.put("S", String.valueOf((int) 2));
		parameters.put("Cn", String.valueOf((int) 255));
		parameters.put("K", String.valueOf((int) 1));

		parameters.put("Bandwidth", "100");
		parameters.put("Mean", "2");
		parameters.put("Std", "1");

		parameters.put("Optical path", "WP");


		ArrayList<String> OPS= new ArrayList<String>();
		OPS.add("WP");
		OPS.add("VWP");
		

		ArrayList<String> virts = new ArrayList<String>();
		//virts.add("VN1");
		//virts.add("VN");
		virts.add("SFC");
		//virts.add("VN4");
		//virts.add("SFC");
		//virts.add("SFC2");
		//virts.add("SFC3");
		
		ArrayList<Boolean> groom = new ArrayList<Boolean>();
		groom.add(true);
		groom.add(false);

		ArrayList<String> groom_string = new ArrayList<String>();
		groom_string.add("G");
		groom_string.add("NGR");

		//for (int o =0;  o<OPS.size();o++)
		//for (int g=0; g<groom.size();g++)
		{
			//parameters.put("grooming", String.valueOf((boolean) groom.get(g)));
			
			for (int v=0; v<virts.size();v++)
			{
				parameters.put("Virtualization", virts.get(v));
				
				
				FileWriter fileWriter = new FileWriter(".\\N7\\"+ "G"+"_"+"VWP" +"_Results_PAL_ILP_"+parameters.get("Virtualization")+"_"+parameters.get("Bandwidth")+"_M"+parameters.get("Mean")+"S"+parameters.get("Std")+".txt",true);
				PrintWriter printWriter = new PrintWriter(fileWriter);
				for (int is=3; is<6;is++)//isolation
					for (int i=1; i<11; i=i+1) //slices
					{
						for (int sim =1; sim< 2; sim++) //simulations
						{ 
							//try {
								Testable_reliableSlicingPAL_NODT Algorithm = new Testable_reliableSlicingPAL_NODT();
								Algorithm.getParameters();

								parameters.put("seed", String.valueOf(sim));	
								parameters.put("isolation", String.valueOf(is));
								parameters.put("S", String.valueOf(i));



								System.out.println(String.format("Testing %1$d slices",i));
								//System.out.println(groom_string.get(g)+" "+OPS.get(o)+" "+virts.get(v)  );
								ReliabilityResult res=	Algorithm.executeAlgorithm(netPlan, parameters, null);

								printWriter.println(res.toString());//write to a file
								//printWriter.println(String.format("%1$d,%2$d,%3$f,%4$f,%5$d,%6$f", res.isolation, res.nslices,res.time,res.wlchann,res.act_nodes,res.load));
								//write to a file


								printWriter.flush();
							//}
							//catch (com.net2plan.interfaces.networkDesign.Net2PlanException e) {System.out.println("!!!--!!!-- NET2PLAN EXCEPTION--!!!--!!!"); sim--;}
							//catch (java.lang.IllegalArgumentException e2) { System.out.println("!!!--!!!-- GENERAL EXCEPTION--!!!--!!!");sim--; }
							//catch ( com.jom.JOMException e3) {System.out.println("!!!--!!!-- JOIM EXCEPTION--!!!--!!!");sim--;}
							//catch (IndexOutOfBoundsException e4) {System.out.println("!!!--!!!-- IOOB EXCEPTION--!!!--!!!");sim--;}
						}
					}
				printWriter.close();
			}
			//Files.write(file, lines, StandardCharsets.UTF_8, StandardOpenOption.APPEND);
			//netPlan.saveToFile(new File("RESULT.n2p"));
		}
	}

}