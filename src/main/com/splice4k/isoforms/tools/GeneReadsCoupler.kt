package com.splice4k.isoforms.tools


import com.splice4k.base.Exons
import com.splice4k.base.Genes
import com.splice4k.base.GenomicLoci
import com.splice4k.index.SJIndex
import com.splice4k.index.AnnotationIndex
import com.splice4k.progressbar.ProgressBar
import java.io.File
import java.io.PrintWriter


/**
 * @author Zhang Yiming
 * @since 2018.09.30
 * @version 20180930
 * 开始对基因和reads进行配对，这个思路比较简单
 */

class GeneReadsCoupler(
    bamIndex: SJIndex,
    reference: AnnotationIndex,
    private val overlapLevel: Double = 90.0
) {
    private val bamIndex = bamIndex.transcripts.sorted()
    private val reference = reference.genes.sorted()
    private val referenceTranscripts = mutableMapOf<String, MutableList<Genes>>()
    private val results = mutableMapOf<Genes, MutableMap<Genes, Int>>()


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

        val pb = ProgressBar(total = this.bamIndex.size.toLong(), message = "testing")
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
            this.getUniqReads(gene = k, reads = v)?.let {

                this.results.putAll(it)
            }
        }
    }


    /**
     * 判断基因和reads是否同源，能够配对
     * @param reads BAM文件中的外显子
     * @param reference 注释文件中的外显子
     * @return 是否有至少一个外显子符合要求，true -> 是，false -> 不是
     */
    private fun isExonMatch( reads: List<Exons>, reference: List<Exons> ): Boolean {
        var i = 0; var j = 0

        while ( i < reads.size && j < reference.size ) {
            val currentReads = reads[i]
            val currentRef = reference[j]

            when {
                currentReads.isUpStream(currentRef) -> j++
                currentReads.isDownStream(currentRef) -> i++
                else -> {
                    if ( currentReads.overlapPercent(currentRef, all = true) > this.overlapLevel ) {
                        return true
                    }
                    j++
                }
            }
        }

        return false
    }


    /**
     * 根据注释，看看能不能将组装的isoforms和注释中的transcripts匹配起来
     *
     */
    private fun matchIsoformsTranscripts( gene: Genes, isoforms: Genes ): Genes {
        val transcripts = this.referenceTranscripts[gene.geneId]!!.sorted()

        var j = 0

        while ( j < transcripts.size ) {
            val currentTran = transcripts[j]

            if (
                    isoforms.overlapPercent(currentTran, all = true) > this.overlapLevel &&
                            this.isExonMatch(isoforms.exons, currentTran.exons)
                    ) {
                isoforms.parent = currentTran.transcriptId
                break
            }

            j++
        }
        return isoforms
    }


    /**
     * 将同基因上的reads融合成特定的独一无二的isoforms，确定为同一个transcripts的融合起来
     * @param gene 注释的基因
     * @param reads 匹配到基因上的所有reads
     * @return 返回一个字典，基因，transcript和count构成的matrix
     */
    private fun getUniqReads( gene: Genes, reads: List<Genes> ): Map<Genes, MutableMap<Genes, Int>>? {
        return when ( reads.size ) {
            0 -> null
            1 -> mapOf(gene to mutableMapOf(reads[0] to 1) )
            else -> {

                val results = mutableMapOf(
                        gene to mutableMapOf<Genes, Int>()
                )

                this.classifyByGenes(reads).forEach {
                    val mergedReads = this.matchIsoformsTranscripts(gene = gene, isoforms = this.mergeGenes(it))

                    val tmpMap = mutableMapOf( mergedReads to it.size)

                    results[gene]!!.putAll(tmpMap)
                }
                return results
            }
        }
    }


    /**
     * 对同一个基因上的所有reads进行分类成不同的isoforms
     * @param reads 某基因上所有的reads
     * @return 分类之后的reads
     */
    private fun classifyByGenes( reads: List<Genes> ): MutableList<MutableList<Genes>> {
        val results =  mutableListOf<MutableList<Genes>>()
        var last = 0
        when(reads.size) {
            1 -> results.add(reads.toMutableList())
            else -> {
                for ( i in 0..(reads.size - 2) ) {
                    if ( !this.isSameLoci(reads[i], reads[i + 1]) ) {
                        this.classifyByExons( reads = reads.subList(last, i + 1), results = results )
                        last = i + 1
                    }
                }

                this.classifyByExons( reads = reads.subList(last, reads.size), results = results )
            }
        }
        return results
    }


    /**
     * 对归属于同一个范围内的reads，再按照exons进行细化分类isoforms
     * @param reads 属于同一个
     */
    private fun classifyByExons( reads: List<Genes>, results: MutableList<MutableList<Genes>> ) {
        when ( reads.size ) {
            1 -> results.add(reads.toMutableList())
            else -> {
                val exons = mutableMapOf<Exons, Int>()
                val tmpExons = mutableListOf<Exons>()
                reads.asSequence().forEach { tmpExons.addAll(it.exons) }
                tmpExons.sort()

                // generate id for all exons
                var id = 1
                for ( i in 0..(tmpExons.size - 2) ) {
                    exons[tmpExons[i]] = id
                    if ( !this.isSameLoci(tmpExons[i], tmpExons[i + 1]) ) {
                        id++
                    }
                }
                exons[tmpExons.last()] = id

                // assign exon id to reads and classification
                val readsId = mutableMapOf<String, MutableList<Genes>>()
                for ( r in reads ) {
                    val exonIdString = r.exons
                            .asSequence()
                            .map { exons[it]!! }
                            .toMutableList()
                            .asSequence()
                            .sorted()
                            .distinct()
                            .joinToString(",")

                    if ( readsId.containsKey(exonIdString) ) {
                        readsId[exonIdString]!!.add(r)
                    } else {
                        readsId[exonIdString] = mutableListOf(r)
                    }
                }

                results.addAll(readsId.values)
            }
        }
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
    private fun mergeGenes( reads: List<Genes> ): Genes {
        val starts = mutableListOf<Int>()
        val ends = mutableListOf<Int>()
        val exons = mutableListOf<Exons>()

        reads.forEach {
            starts.add(it.start)
            ends.add(it.end)
            exons.addAll(it.exons)
        }

        val tmpGenes = Genes(
                chromosome = reads[0].chromosome,
                start = starts.min()!!,
                end = ends.max()!!
        )

        var currentExon: Exons? = null
        for ( i in exons.sorted() ) {
            if ( currentExon == null ) {
                currentExon = i
            }

            if ( this.isSameLoci(currentExon, i) ) {
                currentExon = Exons(
                        chromosome = currentExon.chromosome,
                        start = listOf(currentExon.start, i.start).min()!!,
                        end = listOf(currentExon.end, i.end).max()!!,
                        exonId = ""
                )
            } else {
                tmpGenes.exons.add(currentExon)
                currentExon = i
            }
        }

        tmpGenes.exons.add(currentExon!!)

        return tmpGenes
    }


    /**
     * 将之前保留的特定gtf文件输出出来
     * @param gene 注释基因的范围
     * @param transcripts 匹配到基因上的所有reads构成的不同isoform
     * @return gtf格式的文本
     */
    private fun transcriptToGtf( gene: Genes, transcripts: MutableMap<Genes, Int> ): String {

        val res = mutableListOf(
                "${gene.chromosome}\tSplice4k\tgene\t${gene.start}\t${gene.end}\t.\t${gene.strand}\t.\tgene_id \"${gene.geneId}\"; gene_name \"${gene.geneName}\""
        )
        var index = 0
        for ( (transcript, count) in transcripts ) {
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
                    "cov \"$count\"; " +
                    when (transcript.parent) {
                        "" -> ""
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
                                "exon_id \"$transcript_id.$i\""
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