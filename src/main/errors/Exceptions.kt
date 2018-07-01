package main.errors

class ChromosomeException(override val message: String) : Throwable()

class ExonException(override val message: String) : Throwable()