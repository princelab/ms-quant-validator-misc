#!/usr/bin/env ruby

require 'optparse'
require 'mzml'
require 'gnuplot'

class Array
  def median
    sort[size/2]
  end
end

opt = {
  :type => :median,
  :points => 29,
  :fft_z_start => 0,
  :fft_z_end => -8,
}
opts = OptionParser.new do |op|
  progname = File.basename(__FILE__)
  op.banner = "usage: #{progname} <file>.mzML ..."
  op.separator "output: plots a smoothed chromatogram"
  op.separator ""
  op.separator "options: "
  op.on("-p", "--points <Int>", Integer, "odd number of points (default: #{opt[:points]})") {|v| opt[:points] = v }
  op.on("-s", "--fft-z-start <Int>", Integer, "start of fft, default #{opt[:fft_z_start]}") {|v| opt[:fft_z_start] = v }
  op.on("-e", "--fft-z-end <Int>", Integer, "end of fft, default #{opt[:fft_z_end]}") {|v| opt[:fft_z_end] = v }
  op.on("-t", "--type <String>", "type of smoother: median|fft|none") {|v| opt[:type] = v.to_sym }
  op.on("-o", "--output <String>", "type of image to output: png|svg|gif  etc.") {|v| opt[:ext] = v }
  op.separator ""
  op.separator "example: #{progname} -t fft --ext png file1.mzML file2.mzML"
end
opts.parse!

if ARGV.size == 0
  puts opts
  exit
end

ARGV.each do |file|
  xs = []
  ys = []
  MzML::Doc.open(file) do |mzml|
    mzml.index[:spectrum].keys.each do |id|
      spec = mzml.spectrum(id)
      if mz = spec.mz
        xs << spec.retention_time.to_f
        ys << spec.intensity.first.to_f
      end
    end
  end

  new_xs = []; new_ys = [];
  case opt[:type]
  when :none
    [:fft_z_start, :fft_z_end, :points].each {|v| opt.delete v }
    new_xs = xs
    new_ys = ys
  when :median
    [:fft_z_start, :fft_z_end].each {|v| opt.delete v }
    # inefficient, but simple implementation
    xs.zip(ys).each_cons(opt[:points]) do |ar_of_pairs|
      ys_points = ar_of_pairs.map(&:last)
      xs_points = ar_of_pairs.map(&:first)
      new_xs << xs_points[xs_points.size/2]
      new_ys << ys_points.sort[ys_points.size/2]  # median
    end
  when :fft
    [:points].each {|v| opt.delete v }
    require 'fftw3'
    fft = FFTW3.fft(NArray.to_na(ys))
    start = opt[:fft_z_start]
    fin = opt[:fft_z_end]
    fft[start..fin] = Complex(0)
    normal = FFTW3.ifft(fft)
    new_ys = normal.real.to_a
    new_xs = xs
  end

  Gnuplot.open do |gp|
    Gnuplot::Plot.new(gp) do |plot|
      if ext = opt.delete(:ext)
        plot.terminal ext
        output = [file.chomp(File.extname(file)), opt.map {|pair| pair.join('-')}.join("_")].join("-") + '.' + ext
        plot.output output
      end
      plot.title [file, opt.map {|pair| pair.join(':')}.join(", ")].join(" - ")
      plot.ylabel "ion count"
      plot.xlabel "time (sec)"
      dataset = Gnuplot::DataSet.new( [new_xs, new_ys] ) do |ds|
        ds.with = "lines"
        ds.notitle
      end
      plot.data << dataset
    end
  end
end
