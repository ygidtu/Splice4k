package com.splice4k.isoforms

import com.splice4k.base.Exons
import com.splice4k.base.Reads
import com.splice4k.tools.RelativePos
import java.io.IOException
import java.io.PrintWriter
import java.io.File

/**
 * @author Zhang Yiming
 * @since 20200115
 * @version 20200115
 *
 * 从PacBio or nanopore BAM file to extract reads to unique it
 */



class Uniquor (
        val reads: List<Reads>,
        val n_jobs: Int=1,
        val error: Int=3
) {
    val comp = RelativePos()

    fun unique(): List<Reads> {
        val tmpReads = reads.sorted()

        val res = mutableListOf<Reads>()

        if (tmpReads.isEmpty()) {
            return res
        }

        var tempReads = tmpReads.first()

        for (i in 1.until(tmpReads.size - 1)) {
            val currentReads = tmpReads[i]

            when (comp.isSameReads(tempReads, currentReads))  {
                true -> {

                }
                else -> {
                    res.add(tempReads)
                    tempReads = currentReads
                }
            }
        }

        return res
    }

    private fun exonToString(exon: Exons, transcript_id: String): String {
        return "${exon.chromosome}\tSplice4k\texon\t${exon.start}\t" +
                "${exon.end}\t.\t" +
                "${exon.strand}\t" +
                "transcript_id \"${transcript_id}\"; " +
                "transcript_name \"${transcript_id}\"; " +
                "exon_id \"${transcript_id}.${exon.exonId}\""

    }

    private fun readsToString(reads: Reads): List<String> {
        val res = mutableListOf<String>()

        res.add("${reads.chromosome}\tSplice4k\ttranscript\t${reads.start}\t" +
                "${reads.end}\t.\t" +
                "${reads.strand}\t" +
                "transcript_id \"${reads.name}\"; " +
                "transcript_name \"${reads.name}\"; ")

        reads.getExon().forEach {
            res.add(this.exonToString(it, reads.name))
        }

        return res
    }

    fun writeTo(output: File) {
        try{
            if (!output.absoluteFile.parentFile.exists()) {
                output.absoluteFile.parentFile.mkdirs()
            }

            val writer = PrintWriter(output)

            for (v in this.unique()) {
                this.readsToString(v).forEach {
                    writer.println(it)
                }
            }

            writer.close()
        } catch (err: IOException) {
            println(err.message)
            for (i in err.stackTrace) {
                println(i)
            }
        }
    }
}