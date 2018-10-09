package com.splice4k.sms.tools


import com.splice4k.base.GenomicLoci
import com.splice4k.base.SpliceEvent
import com.splice4k.base.SpliceGraph
import com.splice4k.errors.ChromosomeException
import com.splice4k.index.AnnotationIndex
import com.splice4k.index.SJIndex
import com.splice4k.sms.base.Template
import com.splice4k.tools.CheckAS
import com.splice4k.tools.IdentifyAS
import com.splice4k.tools.PsiOfIR
import org.apache.log4j.Logger
import java.util.concurrent.Executors
import java.util.concurrent.Callable
import java.util.concurrent.Future

/**
 * @since 2018.06.21
 * @version 20180929
 * @author Zhang Yiming
 *
 * 根据基因和Reads的配对情况，找出其中的可变剪接情况
 */

/**
 * @param template 构筑出的模板
 * @param bamIndex 读取的Bam文件的index
 * @param refIndex 读取的reference文件的index
 * @param silent 是否输出详细信息
 * @param overlapOfExonIntron 判断IR时外显子与内含子之间重合程度
 * @param error 判断A3/A5是否存（等）时所需要的阈值
 * @param threads 计算所有的线程数
 */
class SJFinder(
        template: TranscriptsReadsCoupler,
        val bamIndex: SJIndex,
        val refIndex: AnnotationIndex,
        val silent: Boolean,
        val overlapOfExonIntron: Double,
        val error: Int,
        val threads: Int
) {
    private val template = template.templates
    private val logger = Logger.getLogger(SJFinder::class.java.toString())
    val results = hashMapOf<SpliceEvent, MutableList<String>>()
    private val identified = mutableSetOf<String>()

    init {
        this.identifySJ()
    }


    /**
     * 多进程的class，继承自Callable，用以获取完成的分析结果
     * @param gene 某特定基因的id
     * @param template 匹配到该基因上的所有reads构成的模板
     * @param bamIndex 读取的Bam文件的index
     * @param silent 是否输出详细信息
     * @param overlapOfExonIntron 判断IR时外显子与内含子之间重合程度
     * @param error 判断A3/A5是否存（等）时所需要的阈值
     */
    class Run(
            private val gene: String,
            private val template: MutableList<Template>,
            private val bamIndex: SJIndex,
            private val error: Int,
            private val overlapOfExonIntron: Double,
            private val silent: Boolean

    ): Callable<HashMap<SpliceEvent, MutableList<String>>> {
        private val psiOfIR = PsiOfIR()
        private val checkAS = CheckAS()
        private val results = hashMapOf<SpliceEvent, MutableList<String>>()

        /**
         * 找出所有的IR
         * @param template 配对好的模板
         * @return map of SpliceEvent and its corresponding gene, transcript
         */
        private fun findIR(template: Template) {

            val exons = template.template.exons


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

                            tmp.psi = this.psiOfIR.getPsi(
                                    chromosome = template.template.chromosome,
                                    regionStart = exons[i].start - 1,
                                    regionEnd = exons[i].end + 1,
                                    bamFile = this.bamIndex.infile
                            )

                            this.results[tmp] = mutableListOf("${template.template.geneId}\t${template.template.transcriptId}\t${exons[i].exonId}")
                        }
                    }
                }
            }

        }


        /**
         * 重载call function
         * @return 与外部class的结果类别一致
         */
        override fun call(): HashMap<SpliceEvent, MutableList<String>> {
            for ( pair in this.template ) {
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
                    this.checkAS.check(i, pair.template.exons)?.let {
                        this.results[i] = mutableListOf("${this.gene}\t${pair.template.transcriptId}\t${it.exonId}")
                    }
                }

                this.findIR(template = pair)
            }

            return this.results
        }
    }


    /**
     * 识别各种可变剪接类型
     */
    private fun identifySJ() {
        this.logger.info("Finding alternative splicing events")
        // 这个gap就是为了控制输出一个合适的进度条的
        val pool = Executors.newFixedThreadPool(this.threads)
        val futures = mutableListOf<Future<HashMap<SpliceEvent, MutableList<String>>>>()

        for ( (gene, template) in this.template ) {

            val f = pool.submit(Run(
                    gene = gene,
                    template = template,
                    bamIndex = this.bamIndex,
                    error = this.error,
                    overlapOfExonIntron = this.overlapOfExonIntron,
                    silent = this.silent
            ))

            futures.add(f)
        }

        futures.forEach {
            this.results.putAll( it.get() )
        }
        pool.shutdown()


        val identifyAS = IdentifyAS(
                overlapOfExonIntron = this.overlapOfExonIntron,
                bamFile = this.bamIndex.infile
        )

        this.logger.info("Finding alternative splicing events by another algorithm")

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
                    this.results[k] = v.asSequence().map { it.toString() + "\t0" }.toMutableList()
                }
            } else {
                this.results[k] = v.asSequence().map { it.toString() + "\t2" }.toMutableList()
            }
        }

    }
}