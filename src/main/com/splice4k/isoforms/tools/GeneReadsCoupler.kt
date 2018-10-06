package com.splice4k.isoforms.tools


import com.splice4k.base.Exons
import com.splice4k.base.Genes
import com.splice4k.base.GenomicLoci
import com.splice4k.index.AnnotationIndex
import com.splice4k.index.SJIndex
import com.splice4k.isoforms.base.SpliceGraph
import com.splice4k.progressbar.ProgressBar
import java.io.File
import java.io.PrintWriter


/**
 * @author Zhang Yiming
 * @since 2018.09.30
 * @version 20181006
 * 开始对基因和reads进行配对，这个思路比较简单
 */

class GeneReadsCoupler(
    bamIndex: SJIndex,
    reference: AnnotationIndex,
    private val overlapLevel: Double = 90.0
) {
    private val bamIndex = bamIndex.transcripts.sorted()
    private val junctions = bamIndex.data
    private val reference = reference.genes.sorted()
    private val referenceTranscripts = mutableMapOf<String, MutableList<Genes>>()
    private val results = mutableMapOf<Genes, List<Genes>>()


    init {
        reference.transcripts.forEach {
            if (this.referenceTranscripts.containsKey(it.parent)) {
                this.referenceTranscripts[it.parent]!!.add(it)
            } else {
                this.referenceTranscripts[it.parent] = mutableListOf(it)
            }
        }

        this.matchGeneReference()
    }


    /**
     * 将某个特定区域内的reads全部匹配到基因上
     */
    private fun matchGeneReference () {
        val results = mutableMapOf<Genes, MutableList<Genes>>()
        var i = 0; var j = 0

        val pb = ProgressBar(message = "testing")
        while ( i < this.bamIndex.size && j < this.reference.size ) {
            val currentReads = this.bamIndex[i]
            val currentRef = this.reference[j]

            when {
                currentReads.isUpStream(currentRef) -> i++
                currentReads.isDownStream(currentRef) -> { j++; pb.step()}
                else -> {
                    if (
                            currentReads.strand == currentRef.strand &&
                            this.isExonMatch( reads = currentReads.exons, reference = currentRef.exons )
                    ) {
                        if ( results.containsKey(currentRef) ) {
                            results[currentRef]!!.add(currentReads)
                        } else {
                            results[currentRef] = mutableListOf(currentReads)
                        }
                    }
                    i++
                }
            }
        }
        pb.close()

        for ( (k, v) in results ) {
            this.results[k] = this.constructIsoforms(v).map {
                this.matchIsoformsTranscripts(k, it)
            }
        }
    }


    /**
     * 针对每个基因构造图
     * 1. 先收集所有的外显子，整合成交完整的Node
     * 2. 根据junctions，统计weighted
     * 3. 识别所有的isoforms
     * @param genes map到同一个基因的reads
     * @return list of genes，构造好的isoforms
     */
    private fun constructIsoforms( genes: List<Genes> ): List<Genes> {
        val mergedExons = this.mergeExons(genes)

        val junctions = this.junctions["${genes[0].chromosome}${genes[0].strand}"]!!

        val graph = SpliceGraph()

        var i = 0
        while ( i < mergedExons.size - 1 ) {        // for (i, j) in mergedExons.withIndex()
            val juncs = junctions.getSites( mergedExons[i].start to mergedExons[i].end )

            juncs.forEach {
                var j = i + 1

                while ( j < mergedExons.size ) {
                    when  {
                        it.site < mergedExons[j].start -> j ++                        // 位点在exon上游，就下一个位点
                        it.site > mergedExons[j].end -> j = mergedExons.size          // 位点在exon下游，就直接略过了
                        else -> {
                            graph.addEdge( edge = mergedExons[i] to mergedExons[j], weight = it.count )
                            j++
                        }
                    }
                }
            }

            i ++
        }

        val isoforms = graph.getAllIsoforms()

        val res = mutableListOf<Genes>()

        isoforms.forEach {
            val tmpGenes = Genes(
                    chromosome = genes[0].chromosome,
                    start = it.first().start,
                    end = it.last().end
            )

            tmpGenes.exons.addAll( it )

            res.add(tmpGenes)
        }

        return res
    }


    /**
     * 判断基因和reads是否同源，能够配对
     * @param reads BAM文件中的外显子
     * @param reference 注释文件中的外显子
     * @return 是否有至少一个外显子符合要求，true -> 是，false -> 不是
     */
    private fun isExonMatch( reads: List<Exons>, reference: List<Exons> ): Boolean {
        var i = 0; var j = 0
        var matchTimes = 0
        while ( i < reads.size && j < reference.size ) {
            val currentReads = reads[i]
            val currentRef = reference[j]

            when {
                currentReads.isUpStream(currentRef) -> i++
                currentReads.isDownStream(currentRef) -> j++
                else -> {
                    // currentReads.overlapPercent(currentRef, all = true) > this.overlapLevel
                    if (this.isSameLoci( currentReads, currentRef )  ) {
                        currentReads.annotation = currentRef.exonId
                        matchTimes++
                    }
                    i++
                }
            }
        }

        return matchTimes > 0
    }


    /**
     * 根据注释，看看能不能将组装的isoforms和注释中的transcripts匹配起来
     * @param gene 基因
     * @param isoform 组装好的亚型
     * @return Genes 与注释中的内容，做过比对后的基因
     */
    private fun matchIsoformsTranscripts( gene: Genes, isoform: Genes ): Genes {
        val transcripts = this.referenceTranscripts[gene.geneId]!!.sorted()

        var j = 0

        while ( j < transcripts.size ) {
            val currentTran = transcripts[j]

            if (
                    currentTran.overlapPercent(gene) > this.overlapLevel &&
                            this.isExonMatch(isoform.exons, currentTran.exons)
                    ) {
                isoform.parent = currentTran.transcriptId
                break
            }

            j++
        }
        return isoform
    }


    /**
     * 判断这些位点是否为同一个，容许有3bp的误差范围
     * @param first 第一个位点
     * @param second 第二个位点
     * @return boolean true -> 是同一个； false -> 不是同一个
     */
    private fun isSameLoci( first: GenomicLoci, second: GenomicLoci ): Boolean {
        return kotlin.math.abs(first.start - second.start) <= 3 &&
                kotlin.math.abs(first. end - second.end) <= 3
    }


    /**
     * 将聚类到一起的reads融合为一条isofrm
     * @param reads 聚类在一起的reads
     * @return 融合过后的reads
     */
    private fun mergeExons( reads: List<Genes> ): List<Exons> {
        val starts = mutableListOf<Int>()
        val ends = mutableListOf<Int>()
        val exons = mutableListOf<Exons>()

        reads.forEach {
            starts.add(it.start)
            ends.add(it.end)
            exons.addAll(it.exons)
        }

        val results = mutableListOf<Exons>()

        var currentExon: Exons? = null
        for ( i in exons.sorted() ) {
            if ( currentExon == null ) {
                currentExon = i
            }

            if ( currentExon.overlapPercent(i, all = true) >= this.overlapLevel ) {
                currentExon = Exons(
                        chromosome = currentExon.chromosome,
                        start = listOf(currentExon.start, i.start).min()!!,
                        end = listOf(currentExon.end, i.end).max()!!,
                        exonId = ""
                )
            } else {
                results.add(currentExon)
                currentExon = i
            }
        }

        results.add(currentExon!!)

        return results
    }


    /**
     * 将之前保留的特定gtf文件输出出来
     * @param gene 注释基因的范围
     * @param transcripts 匹配到基因上的所有reads构成的不同isoform
     * @return gtf格式的文本
     */
    private fun transcriptToGtf( gene: Genes, transcripts: List<Genes> ): String {

        val res = mutableListOf(
                "${gene.chromosome}\tSplice4k\tgene\t${gene.start}\t${gene.end}\t.\t${gene.strand}\t.\tgene_id \"${gene.geneId}\"; gene_name \"${gene.geneName}\"; "
        )
        var index = 0
        for ( transcript in transcripts ) {
            val transcript_id = "${gene.geneId}.$index"

            res.add(
                    "${gene.chromosome}\t" +
                    "Splice4k\t" +
                    "transcript\t" +
                    "${transcript.start}\t" +
                    "${transcript.end}\t" +
                    ".\t" +
                    "${gene.strand}\t" +
                    ".\t" +
                    "gene_id \"${gene.geneId}\"; " +
                    "gene_name \"${gene.geneName}\"; " +
                    "transcript_id \"$transcript_id\"; " +
                    // "cov \"$count\"; " +
                    when (transcript.parent) {
                        "" -> "ref_transcript \"NA\"; "
                        else -> "ref_transcript \"${transcript.parent}\"; "
                    }
            )

            index ++
            for ( (i, exon) in transcript.exons.withIndex() ) {
                res.add(
                        "${gene.chromosome}\t" +
                                "Splice4k\t" +
                                "exon\t" +
                                "${exon.start}\t" +
                                "${exon.end}\t" +
                                ".\t" +
                                "${gene.strand}\t" +
                                ".\t" +
                                "gene_id \"${gene.geneId}\"; " +
                                "gene_name \"${gene.geneName}\"; " +
                                "transcript_id \"${gene.geneId}.$index\"; " +
                                "exon_id \"$transcript_id.$i\"; " +
                                when (exon.annotation) {
                                    "" -> "ref_exon \"NA\"; "
                                    else -> "ref_exon \"${exon.annotation}\"; "
                                }
                )
            }
        }


        return res.joinToString("\n")
    }


    /**
     * 将结果写出到文件
     * @param output 输出文件
     */
    fun saveTo(output: File) {
        if ( !output.absoluteFile.parentFile.exists() ) {
            output.absoluteFile.parentFile.mkdirs()
        }

        val writer = PrintWriter(output)

        for ( (k, v) in this.results ) {
            writer.println(this.transcriptToGtf(k, v))
        }

        writer.close()

    }
}