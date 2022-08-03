"""Script to download CALIPSO-ICECLOUD from its webpage."""

from esmvaltool.cmorizers.data.downloaders.wget import NASADownloader


def download_dataset(config, dataset, dataset_info, start_date, end_date,
                     overwrite):
    """Download dataset.

    Parameters
    ----------
    config : dict
        ESMValTool's user configuration
    dataset : str
        Name of the dataset
    dataset_info : dict
         Dataset information from the datasets.yml file
    start_date : datetime
        Start of the interval to download
    end_date : datetime
        End of the interval to download
    overwrite : bool
        Overwrite already downloaded files
    """
    downloader = NASADownloader(
        config=config,
        dataset=dataset,
        dataset_info=dataset_info,
        overwrite=overwrite,
    )
    for path in FILES:
        downloader.download_file(path)


SERVER = ("https://asdc.larc.nasa.gov/data/"
          "CALIPSO/LID_L3_Ice_Cloud-Standard-V1-00")

FILES = f"""{SERVER}/2016/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2016-12A.hdf
{SERVER}/2016/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2016-11A.hdf
{SERVER}/2016/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2016-10A.hdf
{SERVER}/2016/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2016-09A.hdf
{SERVER}/2016/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2016-08A.hdf
{SERVER}/2016/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2016-07A.hdf
{SERVER}/2016/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2016-06A.hdf
{SERVER}/2016/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2016-05A.hdf
{SERVER}/2016/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2016-04A.hdf
{SERVER}/2016/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2016-03A.hdf
{SERVER}/2016/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2016-01A.hdf
{SERVER}/2015/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2015-12A.hdf
{SERVER}/2015/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2015-11A.hdf
{SERVER}/2015/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2015-10A.hdf
{SERVER}/2015/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2015-09A.hdf
{SERVER}/2015/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2015-08A.hdf
{SERVER}/2015/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2015-07A.hdf
{SERVER}/2015/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2015-06A.hdf
{SERVER}/2015/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2015-05A.hdf
{SERVER}/2015/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2015-04A.hdf
{SERVER}/2015/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2015-03A.hdf
{SERVER}/2015/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2015-02A.hdf
{SERVER}/2015/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2015-01A.hdf
{SERVER}/2014/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2014-12A.hdf
{SERVER}/2014/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2014-11A.hdf
{SERVER}/2014/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2014-10A.hdf
{SERVER}/2014/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2014-09A.hdf
{SERVER}/2014/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2014-08A.hdf
{SERVER}/2014/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2014-07A.hdf
{SERVER}/2014/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2014-06A.hdf
{SERVER}/2014/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2014-05A.hdf
{SERVER}/2014/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2014-04A.hdf
{SERVER}/2014/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2014-03A.hdf
{SERVER}/2014/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2014-02A.hdf
{SERVER}/2014/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2014-01A.hdf
{SERVER}/2013/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2013-12A.hdf
{SERVER}/2013/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2013-11A.hdf
{SERVER}/2013/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2013-10A.hdf
{SERVER}/2013/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2013-09A.hdf
{SERVER}/2013/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2013-08A.hdf
{SERVER}/2013/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2013-07A.hdf
{SERVER}/2013/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2013-06A.hdf
{SERVER}/2013/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2013-05A.hdf
{SERVER}/2013/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2013-04A.hdf
{SERVER}/2013/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2013-03A.hdf
{SERVER}/2013/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2013-02A.hdf
{SERVER}/2013/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2013-01A.hdf
{SERVER}/2012/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2012-12A.hdf
{SERVER}/2012/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2012-11A.hdf
{SERVER}/2012/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2012-10A.hdf
{SERVER}/2012/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2012-09A.hdf
{SERVER}/2012/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2012-08A.hdf
{SERVER}/2012/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2012-07A.hdf
{SERVER}/2012/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2012-06A.hdf
{SERVER}/2012/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2012-05A.hdf
{SERVER}/2012/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2012-04A.hdf
{SERVER}/2012/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2012-03A.hdf
{SERVER}/2012/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2012-02A.hdf
{SERVER}/2012/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2012-01A.hdf
{SERVER}/2011/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2011-12A.hdf
{SERVER}/2011/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2011-11A.hdf
{SERVER}/2011/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2011-10A.hdf
{SERVER}/2011/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2011-09A.hdf
{SERVER}/2011/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2011-08A.hdf
{SERVER}/2011/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2011-07A.hdf
{SERVER}/2011/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2011-06A.hdf
{SERVER}/2011/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2011-05A.hdf
{SERVER}/2011/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2011-04A.hdf
{SERVER}/2011/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2011-03A.hdf
{SERVER}/2011/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2011-02A.hdf
{SERVER}/2011/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2011-01A.hdf
{SERVER}/2010/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2010-12A.hdf
{SERVER}/2010/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2010-11A.hdf
{SERVER}/2010/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2010-10A.hdf
{SERVER}/2010/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2010-09A.hdf
{SERVER}/2010/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2010-08A.hdf
{SERVER}/2010/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2010-07A.hdf
{SERVER}/2010/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2010-06A.hdf
{SERVER}/2010/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2010-05A.hdf
{SERVER}/2010/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2010-04A.hdf
{SERVER}/2010/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2010-03A.hdf
{SERVER}/2010/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2010-02A.hdf
{SERVER}/2010/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2010-01A.hdf
{SERVER}/2009/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2009-12A.hdf
{SERVER}/2009/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2009-11A.hdf
{SERVER}/2009/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2009-10A.hdf
{SERVER}/2009/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2009-09A.hdf
{SERVER}/2009/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2009-08A.hdf
{SERVER}/2009/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2009-07A.hdf
{SERVER}/2009/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2009-06A.hdf
{SERVER}/2009/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2009-05A.hdf
{SERVER}/2009/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2009-04A.hdf
{SERVER}/2009/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2009-03A.hdf
{SERVER}/2009/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2009-02A.hdf
{SERVER}/2009/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2009-01A.hdf
{SERVER}/2008/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2008-12A.hdf
{SERVER}/2008/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2008-11A.hdf
{SERVER}/2008/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2008-10A.hdf
{SERVER}/2008/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2008-09A.hdf
{SERVER}/2008/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2008-08A.hdf
{SERVER}/2008/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2008-07A.hdf
{SERVER}/2008/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2008-06A.hdf
{SERVER}/2008/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2008-05A.hdf
{SERVER}/2008/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2008-04A.hdf
{SERVER}/2008/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2008-03A.hdf
{SERVER}/2008/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2008-02A.hdf
{SERVER}/2008/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2008-01A.hdf
{SERVER}/2007/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2007-12A.hdf
{SERVER}/2007/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2007-11A.hdf
{SERVER}/2007/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2007-10A.hdf
{SERVER}/2007/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2007-09A.hdf
{SERVER}/2007/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2007-08A.hdf
{SERVER}/2007/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2007-07A.hdf
{SERVER}/2007/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2007-06A.hdf
{SERVER}/2007/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2007-05A.hdf
{SERVER}/2007/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2007-04A.hdf
{SERVER}/2007/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2007-03A.hdf
{SERVER}/2007/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2007-02A.hdf
{SERVER}/2007/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2007-01A.hdf
{SERVER}/2006/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2006-12A.hdf
{SERVER}/2006/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2006-11A.hdf
{SERVER}/2006/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2006-10A.hdf
{SERVER}/2006/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2006-09A.hdf
{SERVER}/2006/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2006-08A.hdf
{SERVER}/2006/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2006-07A.hdf
{SERVER}/2006/CAL_LID_L3_Ice_Cloud-Standard-V1-00.2006-06A.hdf""".split('\n')
