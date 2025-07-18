import typer
from .codon_optimise import codon_optimise

app = typer.Typer()

@app.command()
def main(
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

if __name__ == "__main__":
    app()