package com.splice4k.index

import com.splice4k.base.Exons
import com.splice4k.base.Genes
import com.splice4k.tools.FileValidator
import me.tongfei.progressbar.ProgressBar
import me.tongfei.progressbar.ProgressBarBuilder
import me.tongfei.progressbar.ProgressBarStyle
import org.apache.log4j.Logger
import java.io.File
import java.io.FileInputStream
import java.io.IOException
import java.util.*
import kotlin.system.exitProcess


/**
 * 注释文件index的基础类
 * @version 20181006
 * @author Zhang Yiming
 * @since ???
 */


/**
 * @param infile 输入文件
 * @param smrt 是否为SMRT算法运行
 */
class AnnotationIndex(
        private val infile: File,
        private val smrt: Boolean = false,
        private val iso: Boolean = false
) {
    private val logger = Logger.getLogger(AnnotationIndex::class.java)
    private val fileFormat = FileValidator().check(this.infile)

    // exons(values) separated by chromosome and strand (key)
    val data = mutableMapOf<String, MutableList<Exons>>()

    // list of genes
    val genes = mutableListOf<Genes>()

    // list of transcripts
    val transcripts: MutableList<Genes> = mutableListOf()

    init {
        this.readFromAnnotation()
    }


    private fun getSourceFromGff(info: List<String>): Map<String, String> {
        val results: MutableMap<String, String> = mutableMapOf()

        for (inf in info[0].split(";")) {
            if ( inf == "" ) {
                continue
            }

            val tmp = inf.split("=")

            if ( ":" in tmp[1] ) {
                results[tmp[0]] = tmp[1].split(":")[1]
            } else {
                results[tmp[0]] = tmp[1]
            }
        }

        return results
    }


    private fun getSourceFromGtf(info: List<String>): Map<String, String> {
        val results: MutableMap<String, String> = mutableMapOf()

        for (i in 0..(info.size - 1) step 2) {
            try{
                results[info[i]] = info[i+1]     // get rid of useless characters
                        .replace("\"", "")
                        .replace(";", "")
            } catch ( e: IndexOutOfBoundsException ) {
                continue
            }


        }
        return results
    }


    /**
     * 从Gff文件中读取信息
     */
    private fun readFromAnnotation() {

        try {
            this.logger.info("Reading from ${this.infile}")
            val pbb = ProgressBarBuilder().setStyle(ProgressBarStyle.ASCII).setTaskName("Reading")
            val reader = Scanner(ProgressBar.wrap(FileInputStream(this.infile), pbb))

            // <基因id, 基因本体> 收集基因与它对应的外显子的
            val tmpGenes = mutableMapOf<String, Genes>()
            val geneTranscript = mutableMapOf<String, String>()

            val exons: MutableMap<String, Exons> = mutableMapOf()
            var exonNumber = 1

            // collect transcripts with corresponding exons, to set exon index for ALE, AFE
            val transcriptExon = mutableMapOf<String, MutableList<String>>()

            reader.use {
                while (reader.hasNext()) {
                    val line = reader.nextLine()

                    if ( line.startsWith("#") ) {
                        continue
                    }

                    val lines = line.replace("\n", "").split("\\s+".toRegex())

                    val sources = when(this.fileFormat) {
                        "gtf" ->  this.getSourceFromGtf(lines.subList(8, lines.size))

                        "gff" -> this.getSourceFromGff(lines.subList(8, lines.size))

                        else -> {
                            logger.info("Please check reference file ${this.infile} format, it should be gtf|gff")
                            exitProcess(2)
                        }
                    }


                    if ( this.iso && lines[2].matches(".*gene.*".toRegex(RegexOption.IGNORE_CASE)) ) {
                        val gene = Genes(
                                chromosome = lines[0],
                                start = lines[3].toInt(),
                                end = lines[4].toInt(),
                                strand = lines[6].toCharArray()[0],
                                information = sources
                        )

                        tmpGenes[gene.geneId] = gene

                    }


                    if (lines[2].matches(".*exon.*".toRegex(RegexOption.IGNORE_CASE))) {
                        var exonId = sources["exon_id"] ?: sources["ID"]
                        if ( exonId == null ) {
                            exonId = "${sources["Parent"] ?: sources["transcript_id"]}.$exonNumber"
                            exonNumber ++
                        }

                        val tmp = when( exons.containsKey(exonId) ) {
                            true -> exons[exonId]!!
                            false -> Exons(
                                        chromosome = lines[0],
                                        start = lines[3].toInt(),
                                        end = lines[4].toInt(),
                                        strand = lines[6].toCharArray()[0],
                                        exonId = exonId
                                )
                        }

                        val transcriptId = sources["Parent"] ?: sources["transcript_id"] ?: sources["ID"]!!
                        tmp.source["transcript"]!!.add(transcriptId)
                        tmp.source["gene"]!!.add(geneTranscript[sources["Parent"]] ?: sources["gene_id"] ?: sources["GeneID"] ?: sources["ID"]!!)

                        exons[exonId] = tmp

                        // keep the transcript -> exons relationships
                        val tmpExons = transcriptExon[transcriptId] ?: mutableListOf()
                        tmpExons.add(exonId)
                        transcriptExon[transcriptId] = tmpExons


                        if ( this.iso ) {
                            tmpGenes[sources["gene_id"] ?: sources["GeneID"] ?: geneTranscript[sources["Parent"]!!]  ]!!.exons.add(tmp)
                        }
                    } else if (  // 两种标准，严防不标准的gff文件
                            ("transcript_id" in sources.keys && "exon_id" !in sources.keys ) ||
                            (lines[2].matches(".*(rna|transcript).*".toRegex(RegexOption.IGNORE_CASE)) &&
                                    !lines[2].matches(".*gene.*".toRegex(RegexOption.IGNORE_CASE)))
                    ) {

                        try {
                            geneTranscript[sources["ID"] ?: sources["transcript_id"]!! ] = sources["Parent"] ?: sources["gene_id"] ?: sources["GeneID"]!!
                        } catch (e: NullPointerException) {
                            continue
                        }

                        this.transcripts.add(Genes(
                                chromosome = lines[0],
                                start = lines[3].toInt(),
                                end = lines[4].toInt(),
                                strand = lines[6].toCharArray()[0],
                                information = sources
                        ))

                        exonNumber = 1
                    }
                }

            }

            for (i in this.transcripts) {

                // set exons to transcript
                val exonIds = transcriptExon[i.transcriptId] ?: mutableListOf()
                var tmpExons = mutableListOf<Exons>()
                for ( exonId in exonIds ) {
                    tmpExons.add(exons[exonId]!!)
                }

                i.exons = tmpExons

                if ( tmpExons.size != i.exons.size ) {
                    exitProcess(0)
                }

                // set chrom_pair <-> exon pair to data
                for ( exon in i.exons ) {
                    val chromStrand = "${exon.chromosome}${exon.strand}"
                    tmpExons = this.data[chromStrand] ?: mutableListOf()
                    tmpExons.add(exon)
                    this.data[chromStrand] = tmpExons
                }
            }

            this.genes.addAll(tmpGenes.values)

        }catch (e: IOException) {
            this.logger.error(e.toString())
        }
    }
}