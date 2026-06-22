import io
import sys
import types
import unittest
from unittest import mock


if "pysam" not in sys.modules:
    sys.modules["pysam"] = types.ModuleType("pysam")

from gottcha import gottcha2


class TestGottcha2Cli(unittest.TestCase):
    def test_cli_dispatches_profile(self):
        with mock.patch.object(gottcha2.profile, "main") as profile_main:
            with mock.patch.object(sys, "argv", ["gottcha2", "profile", "-i", "reads.fq"]):
                gottcha2.cli()
        profile_main.assert_called_once_with(["profile", "-i", "reads.fq"])

    def test_cli_prints_version(self):
        buf = io.StringIO()
        with mock.patch.object(sys, "stdout", buf):
            with mock.patch.object(sys, "argv", ["gottcha2", "version"]):
                gottcha2.cli()
        self.assertIn(gottcha2.__version__, buf.getvalue())

    def test_cli_invalid_command_exits(self):
        with mock.patch.object(sys, "argv", ["gottcha2", "badcmd"]):
            with self.assertRaises(SystemExit):
                gottcha2.cli()


if __name__ == "__main__":
    unittest.main()
