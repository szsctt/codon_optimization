import typer
from .codon_optimise import codon_optimise
from .check import check_optimised_sequences
from .gc import check_gc_content

app = typer.Typer()

@app.command()
def opt(
    # Add your command arguments here
    input_file: str = typer.Argument(..., help="Input sequence file"),
    output_file: str = typer.Option(None, help="Output file (default: stdout)"),
    type: str = typer.Option('dna', help="Type of sequence: 'dna' or 'protein'"),
    before: str = typer.Option('', help="Sequence to add before the codon optimised sequence"),
    after: str = typer.Option('', help="Sequence to add after the codon optimised sequence"),
    avoid: list[str] = typer.Option(None, help="List of sequences to avoid in the codon optimised sequence. "
                                               "Use 'N' for any nucleotide, e.g. 'GGCC' or 'GG.NN'. "
                                               "If multiple sequences are provided, they will be avoided in the order "
                                               "they are provided. If a sequence is still present after all codons "
                                               "have been changed, an error will be raised.")
    
):
    """Optimize codons for a given sequence."""
    codon_optimise(input_file, type, before, after, avoid, output_file)

@app.command()
def check(
    reference_file: str = typer.Argument(..., help="Reference fasta file containing original sequences"),
    optimised_files: list[str] = typer.Argument(..., help="One or more optimised fasta files to check"),
):
    """
    Check optimised sequences against reference.
    Compares one or more optimised sequence files against a reference file.
    """
    check_optimised_sequences(reference_file, optimised_files)

@app.command()
def gc(
    seq: str = typer.Argument(..., help="Fasta file containing sequence to check GC content"),
):
    """
    Calculate GC content of a sequence.
    """
    gc = check_gc_content(seq)
    typer.echo(f"GC content: {gc:.2f}%")

if __name__ == "__main__":
    app()