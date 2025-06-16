import sys

from test_utils import MonofonicRun, plot_powerspectra


def main():
    if len(sys.argv) < 3:
        raise RuntimeError("python cli.py plot_powerspectra configs/default_config.conf")
    command = sys.argv[1]
    config = sys.argv[2]
    if "/" in config:
        config = config.split("/")[-1]
    run = MonofonicRun(config)
    if command == "save_reference":
        run.save_as_reference()
    elif command == "plot_powerspectra":
        plot_powerspectra(run)


if __name__ == "__main__":
    main()
