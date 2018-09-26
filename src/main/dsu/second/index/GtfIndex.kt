package dsu.second.index

import dsu.carrier.Exons
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

class GtfIndex(infile: String) : AnnotationIndex (infile) {

    override fun getSource(info: List<String>): Map<String, String> {
        val results: MutableMap<String, String> = mutableMapOf()

        val msgWeNeed = arrayOf(
                "gene_name",
                "gene_id",
                "transcript_name",
                "transcript_id",
                "exon_id",
                "exon_name"
        )

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
                    val sources = this.getSource(lines.subList(8, lines.size))

                    val tmpExon = Exons(
                            chromosome = lines[0],
                            start = lines[3].toInt(),
                            end = lines[4].toInt(),
                            strand = lines[6].toCharArray()[0],
                            exonId = sources["exon_id"]!!
                    )

                    tmpExon.source["transcript"] = sources["transcript_id"]!!
                    tmpExon.source["gene"] = sources["gene_id"]!!

                    val tmp = mutableListOf(tmpExon)
                    val key = "${lines[0]}${lines[6]}"
                    if ( this.data.containsKey(key) ) {
                        this.data[key]!!.addAll(tmp)
                    } else {
                        this.data[key] = tmp
                    }
                }
            }

            reader.close()
            pb.close()
        }catch (e: IOException) {
            this.logger.error(e.toString())
        }

    }
}