from downsample_illumina_pe_reads.downsample_illumina_pe_reads import DownsampleIlluminaPEReads
from unittest.mock import Mock, patch, call
from unittest import TestCase

class TestStringMethods(TestCase):
    def setUp(self) -> None:
        self.dummy_downsampler = DownsampleIlluminaPEReads(None, None, None, None, None)
        self.record_mock_of_sequence_with_size_2 = Mock(sequence="AG")
        self.record_mock_of_sequence_with_size_3 = Mock(sequence="AGC")
        self.record_mock_of_sequence_with_size_4 = Mock(sequence="AGCT")
        self.record_mock_of_sequence_with_size_5 = Mock(sequence="AGCTA")
        self.record_mock_of_sequence_with_size_6 = Mock(sequence="AGCTAC")
        self.read_pair_id_to_number_of_bases = [10, 2, 5, 8]
        self.random_order_of_reads = [1, 2, 0, 3]
        self.reads1_fastx_file = ["read_1", "read_2"]
        self.reads2_fastx_file = ["read_3", "read_4"]
        self.out_reads1_file_write_mock = Mock()
        self.out_reads1_file = Mock(write=self.out_reads1_file_write_mock)
        self.out_reads2_file_write_mock = Mock()
        self.out_reads2_file = Mock(write=self.out_reads2_file_write_mock)

    def test___get_read_pair_id_to_number_of_bases_core___empty_reads___returns_empty_list(self):
        actual = self.dummy_downsampler._get_read_pair_id_to_number_of_bases_core([], [])
        expected = []

        self.assertListEqual(actual, expected)

    def test___get_read_pair_id_to_number_of_bases_core___one_read(self):
        actual = self.dummy_downsampler._get_read_pair_id_to_number_of_bases_core(
            [self.record_mock_of_sequence_with_size_2],
            [self.record_mock_of_sequence_with_size_6])
        expected = [8]

        self.assertListEqual(actual, expected)


    def test___get_read_pair_id_to_number_of_bases_core___two_reads(self):
        actual = self.dummy_downsampler._get_read_pair_id_to_number_of_bases_core(
            [self.record_mock_of_sequence_with_size_6, self.record_mock_of_sequence_with_size_5],
            [self.record_mock_of_sequence_with_size_6, self.record_mock_of_sequence_with_size_4])
        expected = [12, 9]

        self.assertListEqual(actual, expected)


    def test___get_read_pair_id_to_number_of_bases_core___uneven_number_of_reads___smallest_has_no_reads(self):
        actual = self.dummy_downsampler._get_read_pair_id_to_number_of_bases_core(
            [],
            [self.record_mock_of_sequence_with_size_6, self.record_mock_of_sequence_with_size_4])
        expected = []

        self.assertListEqual(actual, expected)

    def test___get_read_pair_id_to_number_of_bases_core___uneven_number_of_reads___smallest_has_one_read(self):
        actual = self.dummy_downsampler._get_read_pair_id_to_number_of_bases_core(
            [self.record_mock_of_sequence_with_size_2],
            [self.record_mock_of_sequence_with_size_6, self.record_mock_of_sequence_with_size_4, self.record_mock_of_sequence_with_size_2])
        expected = [8]

        self.assertListEqual(actual, expected)

    def test___get_read_pair_id_to_number_of_bases_core___uneven_number_of_reads___smallest_has_two_reads(self):
        actual = self.dummy_downsampler._get_read_pair_id_to_number_of_bases_core(
            [self.record_mock_of_sequence_with_size_2, self.record_mock_of_sequence_with_size_3],
            [self.record_mock_of_sequence_with_size_6, self.record_mock_of_sequence_with_size_4, self.record_mock_of_sequence_with_size_2])
        expected = [8, 7]

        self.assertListEqual(actual, expected)


    @patch.object(DownsampleIlluminaPEReads, DownsampleIlluminaPEReads._shuffle_list.__name__)
    def test___get_list_with_random_order_of_read_pair_ids(self, shuffle_list_mock):
        shuffle_list_mock_return = Mock()
        shuffle_list_mock.return_value = shuffle_list_mock_return

        actual = self.dummy_downsampler._get_list_with_random_order_of_read_pair_ids(5)

        shuffle_list_mock.assert_called_once_with([0,1,2,3,4])
        self.assertEqual(shuffle_list_mock_return, actual)

    def test___get_reads_until_bases_are_saturated___number_of_bases_is_zero___read_pair_ids_to_output_is_empty(self):
        self.dummy_downsampler._number_of_bases = 0
        actual = self.dummy_downsampler._get_reads_until_bases_are_saturated(self.read_pair_id_to_number_of_bases,
                                                                             self.random_order_of_reads)
        expected = set()
        self.assertSetEqual(actual, expected)

    def test___get_reads_until_bases_are_saturated___number_of_bases_is_one___read_pair_ids_contains_one_id(self):
        self.dummy_downsampler._number_of_bases = 1
        actual = self.dummy_downsampler._get_reads_until_bases_are_saturated(self.read_pair_id_to_number_of_bases,
                                                                             self.random_order_of_reads)
        expected = {1}
        self.assertSetEqual(actual, expected)

    def test___get_reads_until_bases_are_saturated___number_of_bases_is_two___read_pair_ids_contains_one_id(self):
        self.dummy_downsampler._number_of_bases = 2
        actual = self.dummy_downsampler._get_reads_until_bases_are_saturated(self.read_pair_id_to_number_of_bases,
                                                                             self.random_order_of_reads)
        expected = {1}
        self.assertSetEqual(actual, expected)

    def test___get_reads_until_bases_are_saturated___number_of_bases_is_three___read_pair_ids_contains_two_ids(self):
        self.dummy_downsampler._number_of_bases = 3
        actual = self.dummy_downsampler._get_reads_until_bases_are_saturated(self.read_pair_id_to_number_of_bases,
                                                                             self.random_order_of_reads)
        expected = {1, 2}
        self.assertSetEqual(actual, expected)

    def test___get_reads_until_bases_are_saturated___number_of_bases_is_seven___read_pair_ids_contains_two_ids(self):
        self.dummy_downsampler._number_of_bases = 7
        actual = self.dummy_downsampler._get_reads_until_bases_are_saturated(self.read_pair_id_to_number_of_bases,
                                                                             self.random_order_of_reads)
        expected = {1, 2}
        self.assertSetEqual(actual, expected)

    def test___get_reads_until_bases_are_saturated___number_of_bases_is_eight___read_pair_ids_contains_three_ids(self):
        self.dummy_downsampler._number_of_bases = 8
        actual = self.dummy_downsampler._get_reads_until_bases_are_saturated(self.read_pair_id_to_number_of_bases,
                                                                             self.random_order_of_reads)
        expected = {1, 2, 0}
        self.assertSetEqual(actual, expected)

    def test___get_reads_until_bases_are_saturated___number_of_bases_is_seventeen___read_pair_ids_contains_three_ids(self):
        self.dummy_downsampler._number_of_bases = 17
        actual = self.dummy_downsampler._get_reads_until_bases_are_saturated(self.read_pair_id_to_number_of_bases,
                                                                             self.random_order_of_reads)
        expected = {1, 2, 0}
        self.assertSetEqual(actual, expected)

    def test___get_reads_until_bases_are_saturated___number_of_bases_is_eighteen___read_pair_ids_contains_four_ids(self):
        self.dummy_downsampler._number_of_bases = 18
        actual = self.dummy_downsampler._get_reads_until_bases_are_saturated(self.read_pair_id_to_number_of_bases,
                                                                             self.random_order_of_reads)
        expected = {1, 2, 0, 3}
        self.assertSetEqual(actual, expected)

    def test___get_reads_until_bases_are_saturated___number_of_bases_is_very_high___read_pair_ids_contains_four_ids(self):
        self.dummy_downsampler._number_of_bases = 100000000
        actual = self.dummy_downsampler._get_reads_until_bases_are_saturated(self.read_pair_id_to_number_of_bases,
                                                                             self.random_order_of_reads)
        expected = {1, 2, 0, 3}
        self.assertSetEqual(actual, expected)

    def test___output_reads_core___read_pair_ids_to_output_is_empty___no_reads_are_output(self):
        read_pair_ids_to_output = []

        self.dummy_downsampler._output_reads_core(read_pair_ids_to_output, self.reads1_fastx_file, self.reads2_fastx_file,
                                                  self.out_reads1_file, self.out_reads2_file)

        self.out_reads1_file_write_mock.assert_not_called()
        self.out_reads2_file_write_mock.assert_not_called()


    def test___output_reads_core___read_pair_ids_contains_one_read___one_read_is_output(self):
        read_pair_ids_to_output = [1]

        self.dummy_downsampler._output_reads_core(read_pair_ids_to_output, self.reads1_fastx_file, self.reads2_fastx_file,
                                                  self.out_reads1_file, self.out_reads2_file)

        self.out_reads1_file_write_mock.assert_called_once_with("read_2\n")
        self.out_reads2_file_write_mock.assert_called_once_with("read_4\n")

    def test___output_reads_core___read_pair_ids_contains_two_reads___two_reads_are_output(self):
        read_pair_ids_to_output = [0, 1]

        self.dummy_downsampler._output_reads_core(read_pair_ids_to_output, self.reads1_fastx_file, self.reads2_fastx_file,
                                                  self.out_reads1_file, self.out_reads2_file)

        self.assertEqual(self.out_reads1_file_write_mock.call_count, 2)
        self.out_reads1_file_write_mock.has_calls(call("read_1\n"), call("read_2\n"), any_order=False)
        self.assertEqual(self.out_reads2_file_write_mock.call_count, 2)
        self.out_reads2_file_write_mock.has_calls(call("read_3\n"), call("read_4\n"), any_order=False)

    @patch.object(DownsampleIlluminaPEReads, DownsampleIlluminaPEReads._get_read_pair_id_to_number_of_bases.__name__, return_value=[1,2,3])
    @patch.object(DownsampleIlluminaPEReads, DownsampleIlluminaPEReads._get_list_with_random_order_of_read_pair_ids.__name__)
    @patch.object(DownsampleIlluminaPEReads, DownsampleIlluminaPEReads._get_reads_until_bases_are_saturated.__name__)
    @patch.object(DownsampleIlluminaPEReads, DownsampleIlluminaPEReads._output_reads.__name__)
    def test___downsample_illumina_pe_reads(self, output_reads_mock, get_reads_until_bases_are_saturated_mock,
                                            get_list_with_random_order_of_read_pair_ids_mock,
                                            get_read_pair_id_to_number_of_bases_mock):
        get_list_with_random_order_of_read_pair_ids_return_mock = Mock()
        get_list_with_random_order_of_read_pair_ids_mock.return_value = get_list_with_random_order_of_read_pair_ids_return_mock

        get_reads_until_bases_are_saturated_return_mock = Mock()
        get_reads_until_bases_are_saturated_mock.return_value = get_reads_until_bases_are_saturated_return_mock


        self.dummy_downsampler.downsample_illumina_pe_reads()

        get_read_pair_id_to_number_of_bases_mock.assert_called_once_with()
        get_list_with_random_order_of_read_pair_ids_mock.assert_called_once_with(3)
        get_reads_until_bases_are_saturated_mock.assert_called_once_with([1,2,3],
                                                                         get_list_with_random_order_of_read_pair_ids_return_mock)
        output_reads_mock.assert_called_once_with(get_reads_until_bases_are_saturated_return_mock)