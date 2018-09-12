package dsu.second.index

import java.io.File
import java.util.Scanner

import dsu.carrier.Exons
import dsu.carrier.Genes
import dsu.progressbar.ProgressBar
import java.io.IOException

/**
 * @author Zhangyiming
 * @since 2018.09.03
 * @version 20180903
 * 从Gtf文件中提取所有的exon
 */

class GtfIndex(infile: String) : AnnotationIndex (infile) {

    override fun getSource(info: List<String>): Map<String, String> {
        val results: MutableMap<String, String> = mutableMapOf()

        val msgWeNeed = arrayOf("gene_name", "gene_id", "transcript_name", "transcript_id")

        for (i in 0..(info.size - 1) step 2) {
            if (info[i] !in msgWeNeed) {
                continue
            }
            results[info[i]] = info[i+1]     // get rid of useless characters
                    .replace("\"", "")
                    .replace(";", "")

        }
        return results
    }

    override fun readExon() {

        try {
            val reader = Scanner(File(this.infile))

            val pb = ProgressBar(message = "Reading exon from Gtf")

            while (reader.hasNext()) {
                val line = reader.nextLine()
                val lines = line.split("\\s+".toRegex())
                pb.step()

                if ( line.startsWith("#") ) {
                    continue
                }
                if (lines[2] == "exon") {
                    val tmpExon = Exons(
                            chromosome = lines[0],
                            start = lines[3].toInt(),
                            end = lines[4].toInt(),
                            strand = lines[6].toCharArray()[0],
                            source = this.getSource(lines.subList(8, lines.size))["trancript_id"]!!
                    )
                    this.data.add(tmpExon)
                } else if (lines[2] == "transcript") {
                    val tmpGene = Genes(
                            chromosome = lines[0],
                            start = lines[3].toInt(),
                            end = lines[4].toInt(),
                            strand = lines[6].toCharArray()[0],
                            information = this.getSource(lines.subList(8, lines.size))
                    )

                    this.transcripts[tmpGene.transcriptId] = tmpGene
                }


            }

            reader.close()
        }catch (e: IOException) {
            this.logger.error(e.toString())
        }

    }
}