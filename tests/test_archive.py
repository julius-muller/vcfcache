from pathlib import Path

from vcfstash.utils.archive import extract_cache, file_md5, tar_cache


def test_tar_and_extract_roundtrip(tmp_path):
    cache_dir = tmp_path / "cacheA"
    cache_dir.mkdir()
    (cache_dir / "blueprint").mkdir()
    (cache_dir / "stash").mkdir()
    data_file = cache_dir / "blueprint" / "vcfstash.bcf"
    data_file.write_bytes(b"abc123")

    tar_path = tmp_path / "cacheA.tar.gz"
    tar_cache(cache_dir, tar_path)

    out_dir = tmp_path / "out"
    extracted_root = extract_cache(tar_path, out_dir)

    assert extracted_root.name == "cacheA"
    assert (extracted_root / "blueprint" / "vcfstash.bcf").read_bytes() == b"abc123"


def test_file_md5(tmp_path):
    f = tmp_path / "x.bin"
    f.write_bytes(b"hello")
    assert file_md5(f) == "5d41402abc4b2a76b9719d911017c592"

