"""Isoleucine is reported in scaffolds, having been normalized away for assembly."""

from instanexus.assembly import restore_isoleucine
from instanexus.preprocessing import normalize_sequence, remove_modifications, strip_modifications


class TestStripModifications:
    """`strip_modifications` is `remove_modifications` without the I/L normalization."""

    def test_keeps_isoleucine_where_remove_modifications_collapses_it(self):
        assert remove_modifications("ALIVTQTMK") == "ALLVTQTMK"
        assert strip_modifications("ALIVTQTMK") == "ALIVTQTMK"

    def test_both_strip_the_same_modifications(self):
        for sequence in ("ALKALPM[UNIMOD:35]HIR", "AL(ox)KALPMHIR", "-ALKALPMHIR-"):
            assert normalize_sequence(strip_modifications(sequence)) == remove_modifications(sequence)

    def test_none_survives(self):
        assert strip_modifications(None) is None


class TestRestoreIsoleucine:
    """The reads vote on each I/L position they cover."""

    def test_a_single_read_restores_its_own_isoleucine(self):
        restored, decisions = restore_isoleucine(["AAALVTQTMKAA"], ["ALVTQTMK"], ["AIVTQTMK"])
        assert restored == ["AAAIVTQTMKAA"]
        assert decisions[0]["position"] == 4
        assert decisions[0]["votes_isoleucine"] == 1

    def test_the_majority_decides(self):
        # Three reads of one position, two reading isoleucine.
        restored, _ = restore_isoleucine(["ALVTQTMK"], ["ALVTQTMK"] * 3, ["AIVTQTMK", "AIVTQTMK", "ALVTQTMK"])
        assert restored == ["AIVTQTMK"]

        # The same with the majority the other way.
        restored, _ = restore_isoleucine(["ALVTQTMK"], ["ALVTQTMK"] * 3, ["AIVTQTMK", "ALVTQTMK", "ALVTQTMK"])
        assert restored == ["ALVTQTMK"]

    def test_a_tie_stays_leucine(self):
        # No majority is not a reason to change the residue: leucine is both the
        # status quo and the commoner of the two.
        restored, decisions = restore_isoleucine(["ALVTQTMK"], ["ALVTQTMK"] * 2, ["AIVTQTMK", "ALVTQTMK"])
        assert restored == ["ALVTQTMK"]
        assert decisions == []

    def test_a_position_no_read_covers_stays_leucine(self):
        restored, _ = restore_isoleucine(["ALVTQTMKLLLL"], ["ALVTQTMK"], ["AIVTQTMK"])
        assert restored == ["AIVTQTMKLLLL"]

    def test_positions_are_voted_on_independently(self):
        # Two I/L positions, each decided by the reads covering it alone.
        restored, _ = restore_isoleucine(["ALVTQTMKLR"], ["ALVTQTMKLR"] * 2, ["AIVTQTMKLR", "AIVTQTMKIR"])
        assert restored == ["AIVTQTMKLR"]

    def test_a_read_occurring_twice_votes_at_both_places(self):
        restored, _ = restore_isoleucine(["ALKGGGALK"], ["ALK"], ["AIK"])
        assert restored == ["AIKGGGAIK"]

    def test_a_read_whose_residues_do_not_line_up_cannot_vote(self):
        # A wrong offset is worse than no vote, so a mismatched pair is skipped.
        restored, decisions = restore_isoleucine(["ALVTQTMK"], ["ALVTQTMK"], ["AIVTQTM"])
        assert restored == ["ALVTQTMK"]
        assert decisions == []

    def test_scaffolds_with_no_matching_read_are_returned_unchanged(self):
        restored, decisions = restore_isoleucine(["ALVTQTMK"], ["GGGGG"], ["GGGGG"])
        assert restored == ["ALVTQTMK"]
        assert decisions == []

    def test_nothing_other_than_i_and_l_is_touched(self):
        scaffold = "ALVTQTMKGGWYF"
        restored, _ = restore_isoleucine([scaffold], [scaffold], ["AIVTQTMKGGWYF"])
        assert restored[0].replace("I", "L") == scaffold
