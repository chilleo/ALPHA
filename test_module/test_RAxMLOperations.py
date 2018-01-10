from unittest import TestCase

from module import RAxMLOperations

class TestRAxMLOperations(TestCase):
    def test_taxon_names_getter(self):
        raxOps = RAxMLOperations.RAxMLOperations()
        self.assertEqual(["Seq0", "Seq1", "Seq2", "Seq3"], raxOps.taxon_names_getter("../exampleFiles/2basePhylip.txt"))

    def test_raxml_species_tree(self):
        raxOps = RAxMLOperations.RAxMLOperations()
        print raxOps.raxml_species_tree("../exampleFiles/2basePhylip.txt")

    def test_rooter(self):
        self.fail()

    def test_window_splitter(self):
        self.fail()

    def test_raxml_windows(self):
        self.fail()

    def test_run(self):
        self.fail()
