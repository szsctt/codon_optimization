import pytest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from codon_optimisation.check import check_optimised_sequences

@pytest.fixture
def mock_seqio(monkeypatch):
    """Fixture to mock SeqIO.read with custom sequence records."""
    sequence_records = {}
    
    def _mock_read(filename, format):
        return sequence_records.get(filename)
    
    def _setup_mock(filename, sequence):
        # Create a SeqRecord with an actual Seq object
        seq_obj = Seq(sequence)
        record = SeqRecord(seq_obj, id=f"seq_{filename}")
        sequence_records[filename] = record
    
    monkeypatch.setattr("Bio.SeqIO.read", _mock_read)
    return _setup_mock

def test_matching_sequences(mock_seqio, capsys):
    """Test when all optimized sequences match the reference."""
    # Setup sequences with same translation (GGG, GGT, GGA all translate to G)
    mock_seqio("reference.fasta", "GGG")  # Glycine
    mock_seqio("optimized1.fasta", "GGT")  
    mock_seqio("optimized2.fasta", "GGA")  
    
    check_optimised_sequences("reference.fasta", ["optimized1.fasta", "optimized2.fasta"])
    
    captured = capsys.readouterr()
    assert "All optimised sequences match" in captured.out

def test_non_matching_sequence(mock_seqio):
    """Test when an optimized sequence doesn't match the reference."""
    # Different translations
    mock_seqio("reference.fasta", "ATG")  # Methionine (M)
    mock_seqio("optimized1.fasta", "ATG")  # Methionine (M)
    mock_seqio("optimized2.fasta", "GGT")  # Glycine (G)
    
    with pytest.raises(ValueError) as excinfo:
        check_optimised_sequences("reference.fasta", ["optimized1.fasta", "optimized2.fasta"])
    
    assert "optimized2.fasta" in str(excinfo.value)
    assert "does not match reference" in str(excinfo.value)

def test_single_optimized_file(mock_seqio, capsys):
    """Test with a single optimized file that matches."""
    mock_seqio("reference.fasta", "GGG")
    mock_seqio("optimized.fasta", "GGC")
    
    check_optimised_sequences("reference.fasta", ["optimized.fasta"])
    
    captured = capsys.readouterr()
    assert "All optimised sequences match" in captured.out

def test_empty_optimized_files_list(mock_seqio, capsys):
    """Test with an empty list of optimized files."""
    mock_seqio("reference.fasta", "ATG")
    
    check_optimised_sequences("reference.fasta", [])
    
    captured = capsys.readouterr()
    assert "All optimised sequences match" in captured.out