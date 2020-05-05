####
#
# The MIT License (MIT)
#
# Copyright 2020 Eric Bach <eric.bach@aalto.fi>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is furnished
# to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
####
import abc
from rdkit.Chem import Mol, MolFromSmarts

# TODO: Return RDKit's fingerprint class


class FPrinter(object):
    def __init__(self):
        pass

    @abc.abstractmethod
    def get_binary(self, molecule: Mol) -> list:
        """
        Calculate binary fingerprints for the given molecule.

        :param molecule: RDKit Mol object, molecule for which the fingerprint should be calculated.
        :return: list, fingerprint vector
        """
        pass

    @abc.abstractmethod
    def get_counting(self, molecule: Mol) -> list:
        """
        Calculate counting fingerprints for the given molecule.

        :param molecule: RDKit Mol object, molecule for which the fingerprint should be calculated.
        :return: list, fingerprint vector
        """
        pass


class SubstructureFPrinter(FPrinter):
    def __init__(self, smarts_keys=None):
        self.smarts_keys = self._process_smarts_keys(smarts_keys)

        super(SubstructureFPrinter, self).__init__()

    @staticmethod
    def _process_smarts_keys(smarts_keys: list) -> list:
        if smarts_keys is None:
            smarts_keys = SubstructureFPrinter._load_functional_groups_smarts()
        else:
            assert isinstance(smarts_keys, list)
            assert all([isinstance(key, str) for key in smarts_keys])

        parsed_smarts_keys = SubstructureFPrinter._parse_smarts_keys(smarts_keys)

        return parsed_smarts_keys

    @staticmethod
    def _parse_smarts_keys(smarts_keys: list) -> list:
        olist = []

        for key in smarts_keys:
            olist.append(MolFromSmarts(key))
            if not olist[-1]:
                print("Could not parse SMARTS '%s'." % key)

        return olist


    @staticmethod
    def _load_functional_groups_smarts() -> list:
        smarts_keys = []
        with open("resources/SMARTS_InteLigand.txt", "r") as smarts_file:
            while (line := smarts_file.readline().strip()):
                # Skip comment and empty lines
                if line.startswith("#") or len(line) == 0:
                    continue

                smarts_keys.append(line[line.find(" ") + 1, :])

        return smarts_keys

    def get_counting(self, molecule: Mol) -> list:
        return [len(molecule.GetSubstructMatches(smarts, uniquify=True)) for smarts in self.smarts_keys]

    def get_binary(self, molecule: Mol) -> list:
        return [molecule.HasSubstructMatch(smarts) for smarts in self.smarts_keys]
