/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package commandline;

import java.util.HashMap;
import java.util.Map;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.MutuallyExclusiveGroup;
import net.sourceforge.argparse4j.inf.Namespace;

/**
 *
 * @author jonas
 */
public class ArgumentParser {
    
    
    public static Map<String, Object> parseCommnadLineArguments(String [] args) throws ArgumentParserException{
        
        // build and argument parser and set its properties
        net.sourceforge.argparse4j.inf.ArgumentParser argsParser = ArgumentParsers.newArgumentParser("phenomatch.jar")
                .description("Phenotypic analysis of genes in input regions and "
                        + "phenotypes associated to individuals.")
                .epilog("2014 by Jonas Ibn-Salem <j.ibn-salem@uni-mainz.de>")
                .defaultHelp(true)
                .version("${prog} 0.0.1");    
        //TODO remove hard-coded version. e.g. by this approach:http://stackoverflow.com/questions/2469922/generate-a-version-java-file-in-maven
        
        argsParser.addArgument("-i", "--input-file").required(true)
                .help("input file with genomic regions in TAB separated file format");
        
        argsParser.addArgument("-g", "--genes").required(true).help("Genes in BED like format");
 
        argsParser.addArgument("-O", "--phenotype-ontology").required(true)
               .help("the phenotype ontology in OBO file format");
        argsParser.addArgument("-a", "--annotation-file").required(true)
                .help("phenotype annotation file that maps genes to phenotpye terms");

        argsParser.addArgument("-o", "--output-prefix").required(true)
                .help("prefix of file paths of all output files.");
        // add optional parameters
        argsParser.addArgument("-d", "--domains").required(false)
                .help("Topologically associating domains (TADs) in BED file format. "
                        + "Non-overlapping regions are assumed");        
        argsParser.addArgument("--permut-genes").type(Integer.class).metavar("N")
                .setDefault(0).help("Permute the phenotype annotations of genes "
                        + "associated to human phenotypes N time and run "
                        + "whole analysis as control. ");
        argsParser.addArgument("-v", "--version").action(Arguments.version());
        // build objects to parse the commandline to

        // convert argument into an HashMap
        Namespace ns = null;
        Map<String, Object> argMap = new HashMap<String, Object>(); 
        
        // parse arguments and handle errors
        try{
            ns = argsParser.parseArgs(args);
        }catch (ArgumentParserException e) {
            argsParser.handleError(e);
            System.exit(1);
        }
        
        argMap = ns.getAttrs();

        return argMap;

    }
    
}
