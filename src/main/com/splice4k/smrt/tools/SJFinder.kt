package com.splice4k.smrt.tools


import com.splice4k.base.GenomicLoci
import com.splice4k.base.SpliceEvent
import com.splice4k.base.SpliceGraph
import com.splice4k.errors.ChromosomeException
import com.splice4k.index.AnnotationIndex
import com.splice4k.index.BamIndex
import com.splice4k.smrt.base.Template
import com.splice4k.tools.IdentifyAS
import org.apache.log4j.Logger
import java.io.File
import java.io.PrintWriter
import java.util.concurrent.Executors


/**
 * @since 2018.06.21
 * @version 20180927
 * @author Zhang Yiming
 *
 * 根据基因和Reads的配对情况，找出其中的可变剪接情况
 */


class SJFinder(
        template: GeneReadsCoupler,
        val bamIndex: BamIndex?=null,
        val refIndex: AnnotationIndex,
        val silent: Boolean,
        val overlapOfExonIntron: Double,
        val error: Int,
        val threads: Int
) {
    private val template = template.templates
    private val logger = Logger.getLogger(SJFinder::class.java.toString())
    private val results = hashMapOf<SpliceEvent, MutableList<String>>()
    private val identified = mutableSetOf<String>()


    init {
        this.identifySJ()
    }


    /**
     * 找出所有的IR
     * @param template 配对好的模板
     * @return map of SpliceEvent and its corresponding gene, transcript
     */
    private fun findIR(template: Template) {

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
                        this.results[tmp] = mutableListOf("${template.template.geneId}\t${template.template.transcriptId}\tNA")
                        this.identified.add("${tmp.event}_${tmp.sliceSites}")
                    }
                }
            }
        }

    }


    /**
     * 识别各种可变剪接类型
     */
    private fun identifySJ() {
        this.logger.info("Finding alternative splicing events")
        // 这个gap就是为了控制输出一个合适的进度条的
        val executor = Executors.newFixedThreadPool(this.threads)
        for ( (gene, template) in this.template ) {

            val worker = Runnable {

                for ( pair in template ) {
                    val exons = pair.getReadsExons().iterator()
                    val graph = SpliceGraph(
                            chromosome = pair.template.chromosome,
                            strand = pair.template.strand
                    )

                    for (i in exons) {
                        for (j in 0..(i.size - 1) step 2) {
                            graph.addEdge(start = i[j], end = i[j + 1])
                        }
                    }

                    for (i in graph.identifyAS( error = this.error, silent = this.silent ).iterator()) {
                        this.results[i] = mutableListOf("$gene\t${pair.template.transcriptId}\tNA")
                        this.identified.add("${i.event}_${i.sliceSites}")
                    }

                    this.findIR(template = pair)
                }
            }

            executor.execute(worker)
        }


        executor.shutdown()

        while (!executor.isTerminated) {}

        val identifyAS = IdentifyAS( overlapOfExonIntron = this.overlapOfExonIntron )

        this.logger.info("Finding alternative splicing events by another algorithm")
        bamIndex?.let {
            val tmp = identifyAS.matchEventsWithRef(
                    event = bamIndex.data.values.toList(),
                    annotations = refIndex.data,
                    threads = this.threads,
                    error = this.error,
                    show = false
            )

            for ( (k, v) in  tmp) {
                if ( "${k.event}_${k.sliceSites}" !in this.identified ) {
                    for ( j in v ) {
                        this.results[k] = v.toMutableList()
                    }
                }
            }
        }
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