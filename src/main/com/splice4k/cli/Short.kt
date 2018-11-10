package com.splice4k.cli


import com.github.ajalt.clikt.core.CliktCommand
import com.github.ajalt.clikt.parameters.arguments.argument
import com.github.ajalt.clikt.parameters.arguments.multiple
import com.github.ajalt.clikt.parameters.options.*
import com.github.ajalt.clikt.parameters.types.double
import com.github.ajalt.clikt.parameters.types.file
import com.github.ajalt.clikt.parameters.types.int
import com.splice4k.base.Exons
import com.splice4k.base.SpliceEvent
import com.splice4k.index.AnnotationIndex
import com.splice4k.index.SJIndex
import com.splice4k.tools.CountTable
import com.splice4k.tools.IdentifyAS
import com.splice4k.tools.PSITable
import org.apache.log4j.Logger
import java.io.File
import java.io.PrintWriter
import kotlin.system.exitProcess


/**
 * @author Zhang Yiming
 * @since 2018.09.27
 * @version 20181006
 */


class Short: CliktCommand(help = "Identify alternative splicing events from RNA-seq") {
    private val input by argument(
            help = "Path to input file, multiple files separate by space [BAM|SAM|STAR SJ.out.tab|gmap align|SJ]"
    ).file(exists = true).multiple()

    private val bam by option(
            "-b",
            "--bam",
            help = """
                Path to BAM/SAM files, or the directory contains BAM/SAM files.
                - If input files are BAM/SAM files, this parameter won't work.
                - If auto was set, then this program will auto match BAM/SAM files under the input file directory.
                - If specified path to directory contains BAM/SAM files which name matches the input files, this program will auto match those files
                - If specified BAM/SAM file with this parameter, then this program will calculate PSI of IR using this file
                """
    )

    private val reference by option(
            "-r",
            "--reference",
            help = "Path to reference file"
    ).file(exists = true).required()


    private val output by option(
            "-o",
            "--output",
            help = "Path to output file"
    ).file().required()

    private val outputPSITable by option(
            "--psi-table",
            help = "Output PSI table of same starts and same ends junctions"
    ).flag(default = false)

    private val outputCountTable by option(
            "--count-table",
            help = "Output Count table of Junctions"
    ).flag(default = false)

    private val junctionsFilter by option(
            "-c",
            "--count",
            help = "Filter low abundance junctions [default: 3]"
    ).int().validate {
        it >= 0
    }

    private val overallJunctionFilter by option(
            "--overall-count",
            help="""
                Filter low abundance junctions across all samples.
                Eg: set this parameter to 100, then junctions that total counts in all samples are lower than 100 will filtered.
                Note: this parameter won't disable -c|--count, when this parameter is used, the -c will set to 1 by default
            """
    ).int().validate {
        it >= 0
    }


    private val error by option(
            "-e",
            help = "The error to identify whether AS event exists [default: 3bp]"
    ).int().default(3).validate {
        require( it >= 0 ) {"this value must be positive"}
    }


    private val threads by option(
            "-p",
            "--process",
            help = "Number of processes to use"
    ).int().default(1).validate {
        require( it > 0 && it <= Runtime.getRuntime().availableProcessors())
    }


    private val overlapOfExonIntron by option(
            "--overlap-exon-intron",
            help = "Minimal overlap level between exon with intron required for intron retention identification [default: 90.0]"
    ).double().default(90.0).validate {
        require( it > 0 && it <= 100 ) {"this value must between 0 and 100"}
    }


    private val show by option(
            "-v",
            "-verbose",
            help = "Enable detailed messages"
    ).flag(default = false)


    /**
     * extract matched bam files from input directory
     * @param it: Input files, like SJ.out.tab
     * @param dir: directory store the corresponding BAM/SAM files
     * @return BAM file
     */
    private fun extractBamFile( it: File, dir: File ): File? {
        val pattern = ".*${it.name.replace("[_.]?SJ.out.tab".toRegex(), "")}[._]?(\\w+.)*bam$"
                .toRegex(RegexOption.IGNORE_CASE)

        var res: File? = null

        for ( i in dir.walkTopDown() ) {
            if ( i.name.matches( pattern ) ) {
                res = i.absoluteFile
                break
            }
        }

        return res
    }


    override fun run() {
        val logger = Logger.getLogger(Short::class.java)
        if ( this.output.isDirectory ) {
            logger.error("Please set path of output file [event not exists]")
            exitProcess(0)
        }

        if ( this.input.isEmpty() ) {
            logger.error("input files are required")
        } else {
            val psis = mutableMapOf<SpliceEvent, MutableMap<String, String>>()
            val labels = mutableListOf<String>()
            val results = mutableMapOf<SpliceEvent, MutableList<Exons>>()
            val psiTable = PSITable()

            val junctions = mutableListOf<Pair<SJIndex, File?>>()

            val ref = AnnotationIndex(
                    infile = this.reference.absoluteFile,
                    smrt = false
            )

            val junctionsFilter = when {
                this.junctionsFilter != null -> this.junctionsFilter!!
                this.overallJunctionFilter != null && this.junctionsFilter == null -> 1
                else -> 3
            }

            for ( (idx, it) in this.input.withIndex() ) {
                val sj = SJIndex(
                        infile = it.absoluteFile,
                        filter = junctionsFilter,
                        silent = !this.show,
                        smrt = false
                )


                if ( this.outputPSITable ) {
                    psiTable.addJunctionGraph(sj)
                }


                val bamFile: File?
                if ( sj.fileFormat == "bam" ) {
                    bamFile = it
                } else {
                    bamFile = when ( this.bam ) {
                        null -> this.bam
                        "auto" -> {
                            this.extractBamFile(it, it.absoluteFile.parentFile)
                        }

                        else -> {
                            when {
                                File(this.bam).isDirectory -> {
                                    this.extractBamFile(it, File(this.bam))
                                }
                                File(this.bam).isFile -> File(this.bam)
                                "," in this.bam!! -> {
                                    val files = this.bam!!.split(",")

                                    if ( idx >= files.size ) {
                                        File(files[idx - files.size])
                                    } else {
                                        File(files[idx])
                                    }
                                }
                                else -> null
                            }
                        }
                    }
                }

                junctions.add(Pair(sj, bamFile))

            }


            if ( this.outputCountTable ) {
                val countTable = CountTable()

                countTable.writeTo(
                        this.output,
                        junctions.map { it.first }
                )
            }

            var overallFiltered: HashSet<String>? = null
            if ( this.overallJunctionFilter != null ) {
                val countTable = CountTable()

                overallFiltered = countTable.filter(
                        junctions.map { it.first },
                        this.overallJunctionFilter!!
                )
            }


            // start to calculate AS
            for ( (sj, bamFile) in junctions ) {

                if ( sj.fileFormat == "star" && bamFile != null ) {
                    logger.info( "${sj.infile.name} -> ${bamFile.name}" )
                }

                val identifyAS = IdentifyAS(
                        overlapOfExonIntron = this.overlapOfExonIntron,
                        bamFile = bamFile
                )

                logger.info("Predicting Alternative Splicing events of ${sj.infile.name}")
                val data = identifyAS.matchEventsWithRef(
                        event = sj.getJunctionGraph(overallFiltered),
                        annotations = ref,
                        error = this.error,
                        threads = this.threads,
                        show = this.show
                )

                labels.add( sj.infile.name )
                for ( (k, values) in data ) {
                    if ( results.containsKey(k) ) {
                        results[k]!!.addAll(values)
                    } else {
                        results[k] = values.toMutableList()
                    }

                    if ( psis.containsKey(k) ) {
                        psis[k]!![labels.last()] = k.getPsi()
                    } else {
                        psis[k] = mutableMapOf(labels.last() to k.getPsi())
                    }
                }
            }


            if ( !this.output.absoluteFile.parentFile.exists() ) {
                this.output.absoluteFile.parentFile.mkdirs()
            }

            val writer = PrintWriter(this.output)
            val tmpResults = mutableSetOf<String>()

            writer.println(
                    "#spliceRange\t" +
                    "spliceType\t" +
                    "subtype\t" +
                    "spliceSites\t" +
                    "isNovel\t" +
                    "gene\t" +
                    "transcript\t" +
                    "exon\t" +
                    labels.joinToString("\t")
            )

            for ((key, values) in results ) {
                val psi = mutableListOf<String>()
                for ( label in labels ) {
                    psi.add(psis[key]!![label] ?: "NA")
                }

                val gene = mutableSetOf<String>()
                val transcript = mutableSetOf<String>()
                val exon = mutableSetOf<String>()

                for (v in values) {
                    v.source["gene"]?.let {
                        gene.addAll(it)
                    }

                    v.source["transcript"]?.let {
                        transcript.addAll(it)
                    }

                    exon.add(v.exonId)
                }

                if ( gene.isEmpty() ) gene.add("NA")
                if ( transcript.isEmpty() ) transcript.add("NA")
                if ( exon.isEmpty() ) exon.add("NA")

                tmpResults.add(
                        "$key\t" +
                        "${if (key.isNovel) 1 else 0}\t" +
                        "${gene.joinToString(",")}\t" +
                        "${transcript.joinToString(",")}\t" +
                        "${exon.joinToString(",")}\t" +
                        psi.joinToString("\t")
                )

            }
            writer.print(tmpResults.asSequence().sorted().distinct().joinToString("\n"))
            writer.close()


            if ( this.outputPSITable ) {
                psiTable.writeTo( this.output.absoluteFile )
            }
        }
    }
}