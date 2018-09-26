package dsu.third.template


import dsu.carrier.GenomicLoci
import dsu.carrier.SpliceEvent
import dsu.carrier.SpliceGraph
import dsu.errors.ChromosomeException
import dsu.progressbar.ProgressBar
import dsu.second.identifier.IdentifyAS
import dsu.third.carrier.Template
import dsu.third.extractor.BamExtractor
import dsu.third.extractor.Extractor
import org.apache.log4j.Logger
import java.io.File
import java.io.PrintWriter


/**
 * @since 2018.06.21
 * @version 20180926
 * @author Zhang Yiming
 *
 * 根据基因和Reads的配对情况，找出其中的可变剪接情况
 */


class SJFinder(
        template: GeneReadsCoupler,
        val bamIndex: BamExtractor?=null,
        val refIndex: Extractor,
        val silent: Boolean,
        val overlapOfExonIntron: Double,
        val error: Int
) {
    private val template = template.templates
    private val logger = Logger.getLogger(SJFinder::class.java.toString())
    private val results = hashMapOf<SpliceEvent, MutableList<String>>()
    val identified = mutableSetOf<String>()


    init {
        this.identifySJ()
    }


    /**
     * 找出所有的IR
     * @param template 配对好的模板
     * @return map of SpliceEvent and its corresponding gene, transcript
     */
    private fun findIR(template: Template): Map<SpliceEvent, MutableList<String>> {
        val res = mutableMapOf<SpliceEvent, MutableList<String>>()
        val exonSites = template.template.exons
        val exons = mutableListOf<GenomicLoci>()
        for ( i in 1..(exonSites.size - 2) step 2 ) {
            exons.add(GenomicLoci(
                    start = exonSites[i],
                    end = exonSites[i + 1]
            ))
        }

        val junctions = mutableListOf<GenomicLoci>()

        for ( i in 0..(template.reads.size - 2) ) {
            try {
                junctions.add(GenomicLoci(
                        start = template.reads[i].end,
                        end = template.reads[i + 1].start
                ))
            } catch (e: ChromosomeException) {
                continue
            }

        }

        var i = 0; var j = 0

        while ( i < exons.size && j < junctions.size ) {
            when {
                exons[i].isUpStream(junctions[j]) -> j++
                exons[i].isDownStream(junctions[j]) -> i++
                else -> {
                    if ( exons[i].overlapPercent(junctions[j]) > this.overlapOfExonIntron ) {
                        val tmp = SpliceEvent(
                                chromosome = template.template.chromosome,
                                start = junctions[j].start,
                                end = junctions[j].end,
                                strand = template.template.strand,
                                sliceSites = mutableListOf(
                                        exons[i].end,
                                        exons[i].start,
                                        junctions[j].start,
                                        junctions[j].end
                                ),
                                event = "IR"
                        )
                        res[tmp] = mutableListOf("${template.template.geneId}\t${template.template.transcriptId}\tNA")
                        this.identified.add("${tmp.event}_${tmp.sliceSites}")
                    }
                }
            }
        }

        return res
    }

    /**
     * 识别各种可变剪接类型
     */
    private fun identifySJ() {
        this.logger.info("Finding Alternative Splicing events")
        // 这个gap就是为了控制输出一个合适的进度条的
        val pb = ProgressBar(this.template.size.toLong(), "Splice events identifying at")
        for ( (gene, template) in this.template ) {
            pb.step()
            for ( pair in template ) {
                val graph = SpliceGraph(
                        chromosome = pair.template.chromosome,
                        strand = pair.template.strand
                )

                for (i in pair.getReadsExons() ) {
                    for ( j in 0..(i.size - 1) step 2 ) {
                        graph.addEdge(start = i[j], end = i[j + 1])
                    }
                }

                for ( i in  graph.identifyAS(this.silent)) {
                    this.results[i] = mutableListOf("$gene\t${pair.template.transcriptId}\tNA")
                    this.identified.add("${i.event}_${i.sliceSites}")
                }

                this.results.putAll(this.findIR(template = pair))
            }
        }
        pb.close()

        val identifyAS = IdentifyAS(
                overlapOfExonIntron = this.overlapOfExonIntron,
                distanceError = this.error,
                silent = this.silent
        )

        bamIndex?.let {
            for ( (k, v) in identifyAS.matchEventsWithRef(bamIndex, refIndex) ) {
                if ( "${k.event}_${k.sliceSites}" !in this.identified ) {
                    for ( j in v ) {
                        this.results[k] = v.toMutableList()
                    }
                }
            }
        }

        /*
        if ( bamIndex != null ) {
            this.logger.info("Finding novel AS")

            for (i in bamIndex.graph.values ) {
                for ( j in i.identifyAS(silent) ) {
                    if ( !this.results.containsKey(j) ) {
                        this.results[j] = "NA\tNA"
                    }
                }
            }
        }
        */

    }


    /**
     * 保存至文件
     * @param outfile 输出文件的路径
     * @return
     */
    fun saveTo(outfile: String) {

        if ( !File(outfile).parentFile.exists() ) {
            File(outfile).parentFile.mkdirs()
        }

        val writer = PrintWriter(File(outfile))

        for ( (k, v) in this.results ) {
            for ( j in v) {
                if ( j != "NA\tNA\tNA" ) {
                    writer.println("$k\t$j")
                }
            }
        }

        writer.close()

    }
}