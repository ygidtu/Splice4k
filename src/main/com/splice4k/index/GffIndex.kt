package com.splice4k.index

import com.splice4k.base.Exons
import com.splice4k.base.Genes
import com.splice4k.progressbar.ProgressBar
import java.io.File
import java.io.IOException
import java.lang.NullPointerException
import java.util.*
import kotlin.system.exitProcess


/**
 * @author Zhangyiming
 * @since 2018.09.03
 * @version 20180903
 * 从Gtf文件中提取所有的exon
 */

class GffIndex(
        infile: String,
        smrt: Boolean = true
) : AnnotationIndex ( infile, smrt ) {

    override fun getSource(info: List<String>): Map<String, String> {
        val results: MutableMap<String, String> = mutableMapOf()

        for (inf in info[0].split(";")) {
            val tmp = inf.split("=")

            if ( ":" in tmp[1] ) {
                results[tmp[0]] = tmp[1].split(":")[1]
            } else {
                results[tmp[0]] = tmp[1]
            }
        }

        return results
    }

    override fun readExon() {

        try {
            val reader = Scanner(File(this.infile))

            this.logger.info("Reading from ${this.infile}")
            val pb = ProgressBar(message = "Reading exon from Gff")

            val geneTranscript = mutableMapOf<String, String>()

            val exons: MutableMap<String, MutableList<Int>> = mutableMapOf()


            while (reader.hasNext()) {
                val line = reader.nextLine()

                if ( line.startsWith("#") ) {
                    continue
                }

                val lines = line.replace("\n", "").split("\t".toRegex())
                pb.step()

                val sources = this.getSource(lines.subList(8, lines.size))

                if (  // 两种标准，严防不标准的gff文件
                    "transcript_id" in sources.keys ||
                    (lines[2].matches(".*(rna|transcript).*".toRegex(RegexOption.IGNORE_CASE)) &&
                                !lines[2].matches(".*gene.*".toRegex(RegexOption.IGNORE_CASE)))
                ) {

                    geneTranscript[sources["ID"]!!] = sources["Parent"]!!

                    if ( this.smrt ) {
                        this.transcripts.add(Genes(
                                chromosome = lines[0],
                                start = lines[3].toInt(),
                                end = lines[4].toInt(),
                                strand = lines[6].toCharArray()[0],
                                information = sources
                        ))
                    }

                } else if (lines[2] == "exon") {

                    val tmp = Exons(
                            chromosome = lines[0],
                            start = lines[3].toInt(),
                            end = lines[4].toInt(),
                            strand = lines[6].toCharArray()[0],
                            exonId = sources["exon_id"] ?: sources["ID"]!!
                    )

                    try{
                        tmp.source["transcript"] = sources["Parent"]!!
                        tmp.source["gene"] = geneTranscript[sources["Parent"]]!!
                    } catch (e: NullPointerException) {
                        println(line)
                        println(sources)
                        exitProcess(0)
                    }


                    val tmpExons = mutableListOf(tmp)
                    val key = "${lines[0]}${lines[6]}"
                    if ( this.data.containsKey(key) ) {
                        this.data[key]!!.addAll(tmpExons)
                    } else {
                        this.data[key] = tmpExons
                    }

                    if ( this.smrt ) {
                        val tmpExon = mutableListOf(tmp.start, tmp.end)

                        if (exons.containsKey(sources["Parent"]!!)) {
                            tmpExon.addAll(exons[sources["Parent"]!!]!!)
                        }
                        exons[sources["Parent"]!!] = tmpExon
                    }
                }
            }

            reader.close()
            pb.close()

            if ( this.smrt ) {
                // 为转录本指定外显子
                for (i in this.transcripts) {
                    if (exons.containsKey(i.transcriptId)) {
                        i.exons.addAll(exons[i.transcriptId]!!)
                    }
                }
            }

        }catch (e: IOException) {
            this.logger.error(e.toString())
        }

    }
}