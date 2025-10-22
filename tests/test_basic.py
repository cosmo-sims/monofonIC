from test_utils import MonofonicRun, compare_runs


def test_default_config():
    run = MonofonicRun("default_config.conf")
    run.run()
    compare_runs(run)
