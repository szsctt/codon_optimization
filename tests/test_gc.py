import pytest
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from codon_optimisation.gc import check_gc_content

@pytest.fixture
def mock_sequence(monkeypatch):
    """Fixture to create a mock sequence for testing."""
    def _setup_mock(sequence):
        mock_record = SeqRecord(Seq(sequence), id="test", description="")
        monkeypatch.setattr("Bio.SeqIO.read", lambda *args, **kwargs: mock_record)
    return _setup_mock

def test_gc_content_zero(mock_sequence):
    """Test GC content calculation for a sequence with 0% GC."""
    mock_sequence("AAAATTTT")
    assert check_gc_content("dummy_path.fasta") == 0.0

def test_gc_content_full(mock_sequence):
    """Test GC content calculation for a sequence with 100% GC."""
    mock_sequence("GGGGCCCC")
    assert check_gc_content("dummy_path.fasta") == 100.0

def test_gc_content_half(mock_sequence):
    """Test GC content calculation for a sequence with 50% GC."""
    mock_sequence("ATGC")
    assert check_gc_content("dummy_path.fasta") == 50.0

def test_gc_content_mixed(mock_sequence):
    """Test GC content calculation for a sequence with 75% GC."""
    mock_sequence("ATGCGCGC")
    assert check_gc_content("dummy_path.fasta") == 75.0

def test_gc_content_mixed_case(mock_sequence):
    """Test GC content calculation with mixed case."""
    mock_sequence("ATGCgcgc")
    assert check_gc_content("dummy_path.fasta") == 75.0