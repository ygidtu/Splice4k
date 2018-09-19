package dsu.second.index

import dsu.carrier.Exons
import dsu.carrier.Genes
import dsu.progressbar.ProgressBar
import java.io.File
import java.io.IOException
import java.util.*


/**
 * @author Zhangyiming
 * @since 2018.09.03
 * @version 20180903
 * 从Gtf文件中提取所有的exon
 */

class GffIndex(infile: String) : AnnotationIndex (infile) {

    override fun getSource(info: List<String>): Map<String, String> {
        val results: MutableMap<String, String> = mutableMapOf()

        for (inf in info[0].split(";")) {
            val tmp = inf.split("=")

            if (tmp[0] == "Parent") {
                results["ParentType"] = when {
                    ":" in tmp[1] -> tmp[1].split(":")[0]
                    else -> tmp[1]
                }

                results["Parent"] = when {
                    ":" in tmp[1] -> tmp[1].split(":")[1]
                    else -> tmp[1]
                }
            } else {
                results[tmp[0]] = tmp[1]
            }
        }

        return results
    }

    override fun readExon() {

        try {
            val reader = Scanner(File(this.infile))

            val pb = ProgressBar(message = "Reading exon from Gff")

            while (reader.hasNext()) {
                val line = reader.nextLine()

                if ( line.startsWith("#") ) {
                    continue
                }

                val lines = line.split("\\s+".toRegex())
                pb.step()

                if (lines[2] == "exon") {
                    val tmpExon = Exons(
                            chromosome = lines[0],
                            start = lines[3].toInt(),
                            end = lines[4].toInt(),
                            strand = lines[6].toCharArray()[0],
                            source = this.getSource(lines.subList(8, lines.size))["Parent"]!!
                    )

                    val tmp = mutableListOf(tmpExon)
                    val key = "${lines[0]}${lines[6]}"
                    if ( this.data.containsKey(key) ) {
                        this.data[key]!!.addAll(tmp)
                    } else {
                        this.data[key] = tmp
                    }
                }
//                } else if (lines[2].matches("(.*rna|.*transcript)".toRegex(RegexOption.IGNORE_CASE))) {
//                    val tmpGene = Genes(
//                            chromosome = lines[0],
//                            start = lines[3].toInt(),
//                            end = lines[4].toInt(),
//                            strand = lines[6].toCharArray()[0],
//                            information = this.getSource(lines.subList(8, lines.size))
//                    )
//
//                    this.transcripts[tmpGene.transcriptId] = tmpGene
//                }


            }

            reader.close()
        }catch (e: IOException) {
            this.logger.error(e.toString())
        }

    }
}