package net.maizegenetics.analysis.distance;

import java.awt.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import javax.swing.*;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.genotypecall.AlleleFreqCache;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.plugindef.PluginParameter;
import net.maizegenetics.util.SimpleTableReport;
import net.maizegenetics.util.TableReport;
import org.apache.log4j.Logger;

import net.maizegenetics.analysis.distance.Step1;
import net.maizegenetics.analysis.distance.PreprocessHMP;

/**
 * @author Terry Casstevens
 * @author Zhiwu Zhang
 * @author Peter Bradbury
 *
 */
public class PEVARPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(PEVARPlugin.class);

    // private PluginParameter<PEVAR_METHOD> myMethod = new PluginParameter.Builder<>("method", PEVAR_METHOD.Centered_IBS, PEVAR_METHOD.class)
    //         .guiName("PEVAR method")
    //         .range(PEVAR_METHOD.values())
    //         .description("The Centered_IBS (Endelman - previously Scaled_IBS) method produces a PEVAR matrix that is scaled to give a reasonable estimate of additive "
    //                 + "genetic variance. Uses algorithm http://www.g3journal.org/content/2/11/1405.full.pdf Equation-13. "
    //                 + "The Normalized_IBS (Previously GCTA) uses the algorithm published here: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3014363/pdf/main.pdf.")
    //         .build();

    private PluginParameter<Boolean> myTrueF1s = new PluginParameter.Builder<>("trueF1s", true, Boolean.class)
            .description("Get True F1s").build();

    public PEVARPlugin(){
        super(null, false);
        trueF1s(false);
    }

    public PEVARPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    protected void preProcessParameters(DataSet input) {
        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);

        if ((alignInList == null) || (alignInList.isEmpty())) {
            throw new IllegalArgumentException("PEVARPlugin: Nothing selected. Please select a genotype.");
        }
    }

    @Override
    public DataSet processData(DataSet input) {

        try {
            if(!trueF1s()){
                printResults(input);
                return null;
            }
        

            List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);
            Datum current = alignInList.get(0);
            GenotypeTable alignment = (GenotypeTable) current.getData();
            String name = current.getName();

            List<Datum> summaryTables = new ArrayList<>();

            SimpleTableReport trueF1s = null;
            if(trueF1s()){
                trueF1s = getTrueF1s(alignment); //genos
            } 

            if(trueF1s!=null){
                summaryTables.add(new Datum(name + "_trueF1s", trueF1s, "True F1s of " + name));
            }

            if (summaryTables.isEmpty()){
                return null;
            }

            DataSet output = new DataSet(summaryTables, this);
            fireDataSetReturned(new PluginEvent(output, PEVARPlugin.class));

            return output;
        
        } finally {
            fireProgress(100);
        }
    }

    public static void printResults(DataSet input){
        if(input == null) return;
        List<Datum> alignInList = input.getDataOfType(GenotypeTable.class);
        if(alignInList.isEmpty()) return;
        printResults(alignInList.get(0));
    }

    public static void printResults(Datum current){
        GenotypeTable  alignment = (GenotypeTable) current.getData();
        String name = current.getName();
        printResults(alignment, name);
    }

    public static void printResults(GenotypeTable alignment, String name){
        GenotypeTable genos = alignment;

        long numSites = alignment.numberOfSites();
        long numTaxa = alignment.numberOfTaxa();

        long totalDiploids = numSites * numTaxa;

        System.out.println("Genotype Table Name: " + name);
        System.out.println("Number of Taxa: " + numTaxa);
        System.out.println("Number of Sites: " + numSites);
        System.out.println("Sites x Taxa: " + totalDiploids);

        System.out.println("Chromosomes...");
        Chromosome[] chromosomes = alignment.chromosomes();
        PositionList positions = alignment.positions();

        for (int i = 0; i < chromosomes.length; i++) {
            int[] startEnd = alignment.firstLastSiteOfChromosome(chromosomes[i]);
            System.out.println(chromosomes[i].getName() + ": start site: " + startEnd[0] + " (" + positions.get(startEnd[0]).getPosition() + ") last site: " + startEnd[1] + " (" + positions.get(startEnd[1]).getPosition() + ") total: " + (startEnd[1] - startEnd[0] + 1));
        }
        System.out.println();

    }

    private static int[] implementGetTrueF1s(GenotypeTable genos){
        String fname = genos.taxaName(0).split("-")[0];
        // System.out.println(fname);
        fname=fname+".hmp.txt";
        String[] args = new String[]{fname};
        ArrayList<String> trueF1List = new ArrayList<String>();
        int[] ascore = null;
        PreprocessHMP n = new PreprocessHMP(args);
        if(n.getVal()) {
            Step1 method1 = new Step1(args);
            ascore = method1.getAutomatedScores();
            ArrayList<Integer> mscore = n.getManualScores();

            // ArrayList<String> trueF1List = new ArrayList<Integer>();
            //IDENTIFY TRUE AND FALSE NEGATIVES!!!!!
//          int tp=0, tn=0, fp=0, fn=0;
            
            float count_1 = 0;
            float count_0 = 0;
            System.out.print("\nmanual: ");
            for(int i = 0; i < mscore.size(); i++) System.out.print(mscore.get(i)+" ");
            System.out.print("\nauto:   ");
            for(int i = 0; i < mscore.size(); i++) System.out.print(ascore[i]+" ");
            int tp=0;
            int tn=0;
            int fp=0; 
            int fn=0;
            System.out.println("\n\nConsidered True F1's (upon automation): ");
            for(int i = 0; i < ascore.length; i++) {
                if(mscore.get(i) == 1 && ascore[i] == 1) {
                    tp++;
                    // System.out.println((n.getFilename().split(".hmp.txt")[0])+"-"+(i+1));
                    String nam = ((n.getFilename().split(".hmp.txt")[0])+"-"+(i+1));
                    trueF1List.add(nam);
                }
                else if(mscore.get(i) == 0 && ascore[i] == 0)   tn++;
                else if(mscore.get(i) == 0 && ascore[i] == 1)   fp++;
                else if(mscore.get(i) == 1 && ascore[i] == 0)   fn++;
            }
            System.out.print("\nTP:"+tp+" TN:"+ tn+" FP:"+fp+" FN:"+ fn+"\n\n");
            // System.exit(0);
         //count yung mga number 1
         int count_ng_1 = 0;
         int count_ng_0 = 0;
         for(int j = 0; j < mscore.size(); j++) {
             if(mscore.get(j) == 1) count_ng_1++; 
             else count_ng_0++;
         }
         Double rate = (double) (count_1 / count_ng_1);
         Double rate0 = (double) (count_0 / count_ng_0);
         System.out.println((int) (rate * 100 ) + "% (F1 ONLY)\n"); //IMPORTANT: INDICATE IN LOG FILE THAT ONLY F1'S ARE CONSIDERED. SELFING LINES ARE NOT COVERED BY THE SCOPE OF THE SPECIAL PROBLEM.
        }else {
            if(n.getFilename() != null && !(n.getFilename().equals("whole_snp_seq.hmp.txt"))) 
                System.out.println("!! For file " + n.getFilename() + ": (DATA_ERR) Data anomalies were encountered. No processing was done.\n");
        }

        return ascore;
    } 

    private SimpleTableReport getTrueF1s(GenotypeTable alignment){
        //IMPLEMENT STEP 1 OF PEVAR HERE

        // int[] result = implementGetTrueF1s(alignment);
        // Step1 method = implementGetTrueF1s(alignment);
        int[] result = implementGetTrueF1s(alignment);
        // System.out.println(method.get)
        // int[] result = new int[]{0,0,0,0,0,0,0,0,0,0};
        String[] firstColumnNames = new String[]{"Taxa Name", "Result"};

        long numTaxa = alignment.numberOfTaxa();

        Object[][] data = new Object[(int)numTaxa][firstColumnNames.length];

        for(int c = 0; c < numTaxa; c++){
            data[c][0] = alignment.taxaName(c);
            data[c][1] = result[c];
        }
        // result=null;

        return new SimpleTableReport("True F1 Results", firstColumnNames, data);
    }

    public TableReport[] runPlugin(GenotypeTable genotype) {
        DataSet input = new DataSet(new Datum("Genotype Table", genotype, null), this);
        DataSet dataSet = performFunction(input);
        TableReport[] result = new TableReport[dataSet.getSize()];
        for (int i = 0; i < dataSet.getSize(); i++) {
            result[i] = (TableReport) dataSet.getData(i).getData();
        }
        return result;
    }

    public Boolean trueF1s(){
        return myTrueF1s.value();
    }

    public PEVARPlugin trueF1s(Boolean value){
        myTrueF1s = new PluginParameter<>(myTrueF1s, value);
        return this;
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = PEVARPlugin.class.getResource("/net/maizegenetics/analysis/images/Kin.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getButtonName() {
        return "Generate True F1s";
    }

    @Override
    public String getToolTipText() {
        return "Obtain True F1's from marker data";
    }

}



/*

public void listFilesForFolder(final File folder) {
    for (final File fileEntry : folder.listFiles()) {
        if (fileEntry.isDirectory()) {
            listFilesForFolder(fileEntry);
        } else {
            System.out.println(fileEntry.getName());
        }
    }
}

final File folder = new File("/home/you/Desktop");
listFilesForFolder(folder);

*/